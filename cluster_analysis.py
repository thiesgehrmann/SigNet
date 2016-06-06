#!/bin/env python

from ibidas import *;
import os;
import numpy as np;

  # Image stuff
import matplotlib.pylab as plt;
import seaborn as sns;

import scipy.spatial.distance as spd;
import fastcluster as fc;

import DTcut as dtcut
import python_utils as utils;

###############################################################################

# If we decide to use the original thresholded string network, use this as W_PPI
# Intuitively, the diffused network will not introduce links that are stronger tha
#W_PPI = python_utils.load_network_tsv('analysis/diffusion_input_network.tsv');


  # Truncated STRING
os.sys.argv = [ os.sys.argv[0], 'analysis/diffusion_input_network.tsv', 'analysis/genes_in_network.dat', 'analysis/diffusion_output_scores.0.00000.tsv', 0.00 ];

  # JACCARD STRING
os.sys.argv = [ os.sys.argv[0], 'analysis/jaccard_network.tsv', 'analysis/genes_in_network.dat', 'analysis/diffusion_output_scores.0.00000.tsv', 0.00 ];

if __name__ == '__main__':

  N_file = os.sys.argv[1];        # Network file
  V_file = os.sys.argv[2];        # Vertices in network file (genes in network)
  S_file = os.sys.argv[3];        # Scores file
  beta   = float(os.sys.argv[4]); # Diffusion value

    # Read 
  print "Reading data"
  W_PPI  = utils.load_network_tsv(N_file);
  genes_in_network = Load(V_file)
  scores = Read(S_file, delimiter='\t').Cast(float)()

    # Prepare and perform the clustering of the network data
  print "Generating dendrogram"
  Y = 1 - spd.squareform(W_PPI)
  D = fc.linkage(Y, 'centroid');
  T = dtcut.prepare_tree(D);

    # Initialize the tree
  print "Initializing DT"
  dt_tree = dtcut.DTCUT(T, W_PPI, genes_in_network.protein(), scores);

    # Initialize the cltTest
  print "Defining test"
  pvalue_thresh = 0.05;
  min_set_size  = 11;
  max_set_size  = sys.maxint;
  stat_test     = dtcut.cltTest(dt_tree, pvalue_thresh);

    # Find the significant clusters
  print "Finding clusters"
  S = dt_tree.test_tree(stat_test, min_set_size, max_set_size);

    # Retrieve information from each cluster
  clusters = dt_tree.get_clusters(S);
  info     = dt_tree.get_clusters_info(S);

    # Save the information from each cluster
  print "Saving..."
  C_file = 'analysis/clusters_output.%0.5f.dat'
  if len(info) > 0:
    C_info = Rep([ (i, h, p, m) for (i, (h, p, members)) in enumerate(info) for m in members ]) / ('cluster_id', 'height', 'qvalue', 'member_id')
    Save(C_info, C_file);
  #fi

#fi


