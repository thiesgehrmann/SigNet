#!/bin/env python

from ibidas import *;
from bx.intervals.intersection import IntervalTree
from scipy import stats as sst
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
#W_PPI = python_utils.load_network_tsv('analysis/diffusion_input_network.tsv');

genes_in_network = Load('analysis/genes_in_network.dat');


for beta in np.linspace(0, 0.03, 10):

  N_file = 'analysis/diffusion_output_network.%0.5f.tsv' % beta;
  S_file = 'analysis/diffusion_output_scores.%0.5f.tsv' % beta;

    # Read 
  W_PPI  = python_utils.load_network_tsv(N_file);
  scores = Read(S_file)()

    # Prepare and perform the clustering of the network data
  Y = 1 - spd.squareform(W_PPI)
  D = fc.linkage(Y, 'single');
  T = dtcut.prepare_tree(D);

    # Initialize the tree
  dt_tree = dtcut.DTCUT(T, W_PPI, gene_in_network.protein(), gene_scores);

    # Initialize the cltTest
  stat_test     = dtcut.cltTest(dt_tree);
  pvalue_thresh = 0.05;
  min_set_size  = 10;
  max_set_size  = 400;

    # Find the significant clusters
  S = dt_tree.test_tree(stat_test, pvalue_thresh, min_set_size, max_set_size);

    # Retrieve information from each cluster
  clusters = current_tree.get_clusters(S);
  info     = current_tree.get_clusters_info(S);

    # Save the information from each cluster
  C_file = 'analysis/clusters_output.%0.5f.dat'
  if len(info) > 0:
    C_info = Rep([ (i, h, p, m) for (i, (h, p, members)) in enumerate(info) for m in members ]) / ('cluster_id', 'height', 'qvalue', 'member_id')
    Save(C_info, C_file);
  #fi

    # Make space, just in case. Delete this shit
  del dt_tree;
  del stat_test;
  del T; 

#efor


