import numpy as np;
from scipy.spatial.distance import squareform

###############################################################################

def diffuse_kernel(N, BETA):

  return_first = False

  beta_kernels = [];

  if not(hasattr(BETA, '__len__')):
    BETA = [ BETA ];
    return_first = True
  #fi

  A          = np.diag(np.sum(N, 1))
  L          = N - A;
  Eval, Evec = np.linalg.eig(L);
  EvecT      = np.linalg.inv(Evec);
  for beta in BETA:
    temp       = beta*Eval;
    temp[temp < -2**10] = -2**10;
    K = Evec * np.diag(np.exp(temp)) * EvecT;
    K[K< 0]    = 0;
    beta_kernels.append(K);
  #efor

  return beta_kernels[0] if return_first else beta_kernels;
#edef


###############################################################################

def random_network(n=20, e=20, s=5):

  L     = np.zeros(n*(n-1)/2.0)
  Ri    = np.random.choice(range(len(L)), size=e, replace=False)
  Rv    = np.array([ min(np.abs(np.random.randn()/2),1) for x in Ri])
  L[Ri] = Rv
  N     = squareform(L)
 
  S      = np.zeros(n);
  RSi    = np.random.choice(range(n), size=s, replace=False)
  RSv    = np.array([ min(np.abs(np.random.randn()/2),1) for x in RSi])
  S[RSi] = RSv;

  return N, S;
#edef

###############################################################################


