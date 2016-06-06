from ibidas import *;
import os;
import sys
import numpy as np;

  # Image stuff
import matplotlib.pylab as plt;
import seaborn as sns;

import scipy.cluster.hierarchy as sph
import scipy.spatial.distance as spd;
import fastcluster as fc;

import python_utils as utils;

########### FUCK THIS SHIT!!!

sys.setrecursionlimit(10000)

betas = [ 0.00000, 0.00330, 0.00500, 0.00670, 0.01000, 0.01330, 0.01500, 0.01670, 0.02000, 0.02330, 0.02500, 0.02670, 0.03000 ]

genes_in_network = Load('analysis/genes_in_network.dat')
network_order = Read('analysis/network_order.tsv').Cast(int)();

scores = [ ]

for b in betas:
  S = Read('analysis/diffusion_output_scores.%0.5f.tsv' % b, delimiter='\t');
  scores = scores + [ (b, name, score) for (name, score) in zip(genes_in_network.protein(), S())];
#efor

R  = Rep(scores).Cast(float, str, float)
DF = utils.rep_to_df(R)

Export(R, 'analysis/all_diffused_scores.tsv')

DF.columns = ['beta', 'protein', 'score']

DF = DF.pivot('beta', 'protein', 'score');

colnames = DF.columns.tolist()
new_order = [ colnames[x-1] for x in network_order];

DF_ordered = DF[new_order]

sns.heatmap(DF_ordered)

