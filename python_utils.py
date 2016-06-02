import pandas as pd;
from ibidas import *;
import numpy as np;

###############################################################################

def rep_to_df(rep):
  names = rep.Names;
  data = rep();

  df = pd.DataFrame({ n : d for (n,d) in zip(names, data) })

  return df[names];
#edef

###############################################################################
  # From http://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail

def tail( f, lines=20 ):
    total_lines_wanted = lines

    BLOCK_SIZE = 1024
    f.seek(0, 2)
    block_end_byte = f.tell()
    lines_to_go = total_lines_wanted
    block_number = -1
    blocks = [] # blocks of size BLOCK_SIZE, in reverse order starting
                # from the end of the file
    while lines_to_go > 0 and block_end_byte > 0:
        if (block_end_byte - BLOCK_SIZE > 0):
            # read the last block we haven't yet read
            f.seek(block_number*BLOCK_SIZE, 2)
            blocks.append(f.read(BLOCK_SIZE))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            blocks.append(f.read(block_end_byte))
        lines_found = blocks[-1].count('\n')
        lines_to_go -= lines_found
        block_end_byte -= BLOCK_SIZE
        block_number -= 1
    all_read_text = ''.join(reversed(blocks))
    return '\n'.join(all_read_text.splitlines()[-total_lines_wanted:])

###############################################################################

def save_network_tsv(N, output_file):
  fd = open(output_file, 'w');
  n_nodes = N.shape[0];
  for (i,j) in zip(*np.nonzero(N)):
    fd.write("%d\t%d\t%f\n" % (i+1, j+1, N[i,j]));
  #efor
  fd.write("%d\t%d\t%f\n" % (n_nodes, n_nodes, 0.0));
  fd.close();
#edef

###############################################################################

def load_network_tsv(input_file):
  fd = open(input_file, 'r');

  n = int(tail(fd, 1).split('\t')[0]);

  N = np.zeros([n, n]);
  for line in fd:
    i, j, w = tuple(line.strip().split('\t'))
    N[i,j] = w;
  #efor
  return np.matrix(N);
#edef


###############################################################################
