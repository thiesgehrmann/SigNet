#!/bin/env python

from ibidas import *;
from bx.intervals.intersection import IntervalTree
from scipy import stats as sst
import os;
import numpy as np;

import python_utils as utils;
import matplotlib.pylab as plt;
import seaborn as sns;

###############################################################################

def index_gff3_genes(G, window=0):

  G = G.GroupBy(_.seqname).Sort(_.start);
  G = G.Get(_.seqname, _.name, _.start, _.end, _.strand).Flat();

  chrs = {};
  for (seqname, name, start, end, strand) in zip(*G()):
    if seqname not in chrs:
      chrs[seqname] = IntervalTree();
    #fi

    chrs[seqname].add(start-window, end+window, (name, start, end));
    print seqname, start, end, name;
  #efor
    
  return chrs;
#edef

###############################################################################

def match_mut_to_gene_flattop(G_index, D):

  D_gene_score = [];

  for i,mut in enumerate(zip(*D())):
    for (name, start, end) in G_index[mut[0]].find(int(mut[1]),int(mut[1])):
       if mut[1] >= start and mut[1] <= end:
         score = sst.norm.pdf(0);
       else:
        score = sst.norm.pdf( float(min(np.abs(start - mut[1]), np.abs(end - mut[1]))) / 25000.0)
       #fi
       D_gene_score.append((i,) + tuple(mut) + (name,score))
  #efor
  
  D_match =  Rep(D_gene_score) / (('i',) + tuple(D.Names) + ('gene','score'))
  return D_match;

#edef

###############################################################################
###############################################################################

  # Calculate mutation scores

#D = Read('data/NCBIBuild37-mm9-p53-p19-wt-Barcode.csv');
#D = D.Cast(str, int, str, int, str, str, str, str, str, str, str, str)

if not(os.path.exists('analysis/G_score.dat')):

  DT = Read('data/NCBIBuildM37_mm9_All_Genotypes_Tissue.csv')
  DT = DT.Cast(str, int, str, int, str, str,str,str,str,str,str,str,str,int,int, str, float)
  
  G = Read('data/Mus_musculus.NCBIM37.67.gff');
  
    # Put geneID value into name value
  G = G.To(_.name, Do=_.Get(G.Get(_.attr).Each(lambda x: x.split(';')[[ i for (i,n) in enumerate(x.split(';')) if 'geneID' in n][0]].split('=')[1] if 'geneID' in x else '')) /'name' )
  G = G.Detect();
  
    # Get gene features
  G_genes = G[_.feature == 'transcript'].GroupBy(_.name).Get(_.seqname[0], _.source[0], _.id.Array().Each(lambda x: ';'.join(x)), _.parent[0], _.name, _.feature[0], _.start.Min(), _.end.Max(), _.score[0], _.strand[0], _.frame[0], _.attr[0])
  G_names = G_genes.Get(_.name, _.Get(_.attr).Each(lambda x: x.split(';')[[ i for (i,n) in enumerate(x.split(';')) if 'gene_name' in n][0]].split('=')[1] if 'gene_name' in x else ''))
  
    # Construct an index of all genes with a window of 50kb, as in Sepideh's work
  genome_index_window = index_gff3_genes(G_genes, window=50000)
  
    # Score each mutation based on its vicinity to a gene (give one score per gene)
  #D_score = match_mut_to_gene_flattop(genome_index_window, D);
  DT_score = match_mut_to_gene_flattop(genome_index_window, DT);
  
    # Score each gene based on the number of mutations (add up the mutation scores)
  G_score = DT_score.GroupBy(_.gene).Get(_.gene, _.score.Sum());
  G_score = (G_names | Match(_.name, _.gene, jointype='left', merge_same='equi') | G_score).ReplaceMissing()
  Save(G_score, 'analysis/G_score.dat')
else:
  G_score = Load('analysis/G_score.dat');
#fi

if not(os.path.exists('analysis/P_score.dat')):
  id_transfer = Load('data/10090.protein.aliases.v10.txt')
  P_score = (id_transfer[_.alias.In(G_score.name)] | Match(_.alias, _.name, jointype='left', merge_same='equi') | G_score).Get(_.string_protein_id.Each(lambda x: x.split('.')[1]).Cast(str), _.alias / 'gene_id', _.attr, _.score)
  P_score = P_score.ReplaceMissing()
  P_score = P_score.GroupBy(_.string_protein_id)[_.score == _.score.Max()].Flat()
  P_score = P_score.Unique(_.string_protein_id)
  P_score = P_score.Copy()
  Save(P_score, 'analysis/P_score.dat');
else:
  P_score = Load('analysis/P_score.dat');
#fi

###############################################################################
###############################################################################
  # Normalize the network scores

if not(os.path.exists('analysis/string_norm.dat')):
  string    = Read('data/10090.protein.links.v10.txt', delimiter=' ')
  string    = string.Cast(str, str, int);
  string    = string.To(_.protein1, Do=_.Each(lambda x: x.split('.')[1]));
  string    = string.To(_.protein2, Do=_.Each(lambda x: x.split('.')[1]));
  min_score = string.combined_score.Min()()
  max_score = string.combined_score.Max()()
  string    = string.To(_.combined_score, Do=_.Each(lambda x: float((x - min_score))/float(max_score)))
  string    = string.Cast(str, str, float);
  string    = string.Copy();
  Export(string, 'data/string_norm.tsv')
  string = Read('data/string_norm.tsv');
  string = string.Cast(str, str, float)
  Save(string, 'analysis/string_norm.dat');
else:
  string = Load('analysis/string_norm.dat');
#fi

###############################################################################
###############################################################################
  # Truncate the network at 0.5
string_thresh = 0.5;

string_trunc = string[_.combined_score > string_thresh].Copy();

###############################################################################
###############################################################################
  # Prepare network for diffusion

  # First, let's select only the genes that are actually in our network
  # (Some mutated genes may not be here)
genes_in_network = ((string_trunc.protein1 | Stack | string_trunc.protein2).Unique() / ('protein',)) | Match(_.protein, _.string_protein_id, jointype='left', merge_same='equi')  | P_score
genes_in_network = genes_in_network.ReplaceMissing().Copy();

  # Create an index of all the genes in the network
gene_scores = genes_in_network.score();
gene_names  = genes_in_network.protein();
gene_dict   = dict((gene_id, i) for (i, gene_id) in enumerate(genes_in_network.protein()))

  # Create the empty matrix
W_PPI = np.zeros((len(gene_dict), len(gene_dict)))

  # And add in the edges
for (p1, p2, score) in zip(*string_trunc()):
  W_PPI[gene_dict[p1], gene_dict[p2]] = score;
  W_PPI[gene_dict[p2], gene_dict[p1]] = score;
#efor

 # Save it
utils.save_network_tsv(W_PPI, 'analysis/diffusion_input_network_sparse.tsv')
Export(genes_in_network.score, 'analysis/diffusion_input_scores.tsv', names=False);
Save(genes_in_network, 'analysis/genes_in_network.dat');

###############################################################################

  # Draw some figures

plt.cla(); plt.clf();
g = sns.distplot(gene_scores[gene_scores > 0], bins=200)
plt.xlim([0, np.max(gene_scores)]);
plt.xlabel("Gene Mutation score");
plt.ylabel("Frequency");
#g.axes.set_yscale('log')
plt.savefig("analysis/figures/gene_mutation_score_histogram.svg");

###############################################################################

  # Prepare the clustering network
string_full = (string_trunc | Stack | string_trunc.Get(1,0,2)) / ('protein1', 'protein2', 'score');
string_full = string_full.To(_.protein1, Do=_.TakeFrom(gene_dict)).To(_.protein2, Do=_.TakeFrom(gene_dict)).Copy();
string_full = string_full / ('protein1', 'protein2', 'score');
neighbors   = dict([ (p1, set(N+[p1])) for (p1, N) in zip(*string_full.GroupBy(_.protein1).Without(_.score)())]);

JACC = np.zeros((len(gene_dict), len(gene_dict)))

n = len(gene_dict)
for i in xrange(0, n-1):
  for j in xrange(i+1, n):
    n_i  = neighbors[i];
    n_j  = neighbors[j];
    n_ij = n_i | n_j;
    JACC[i,j] = float(len(n_i & n_j)) / float(len(n_ij));
    JACC[j,i] = JACC[j,i];
  #efor
#efor

utils.save_network_tsv(JACC, 'analysis/jaccard_network.tsv')
