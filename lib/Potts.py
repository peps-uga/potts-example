######################
## Imports
######################

import numpy as np



######################
##  L2 Potts
######################


def L2_Potts(y,gamma):
	n = y.size
	z = np.zeros(n)

	# FindBestPartition
	######################

	B = - gamma*np.ones(n+1)
	p = np.zeros(n,dtype=int)

	for r in range(n):
		B[r+1] = np.inf
		for l in range(r):
			b = B[l] + gamma + dstar_L2(y[l:r+1])
			if b<=B[r+1]:
				B[r+1] = b
				p[r] = l-1

	# SegmenationFromPartition
	######################

	r = n ; l = p[r-1] ;

	while r>0:
    		for t in range(l+1,r):
        		z[t] = mustar_L2(y[l+1:r])
    		r = l+1 ; l = p[r-1] ;

	return z




######################
##  L1 Potts
######################


def L1_Potts(y,gamma):
	n = y.size
	z = np.zeros(n)

	# FindBestPartition
	######################

	B = - gamma*np.ones(n+1)
	p = np.zeros(n,dtype=int)

	for r in range(n):
		B[r+1] = np.inf
		for l in range(r):
			b = B[l] + gamma + dstar_L1(y[l:r+1])
			if b<=B[r+1]:
				B[r+1] = b
				p[r] = l-1

	# SegmenationFromPartition
	######################

	r = n ; l = p[r-1] ;

	while r>0:
    		for t in range(l+1,r):
        		z[t] = mustar_L1(y[l+1:r])
    		r = l+1 ; l = p[r-1] ;

	return z







######################
## Distance Functions
######################

def dstar_L2(v):
	return np.linalg.norm(v-np.mean(v),2)**2

def mustar_L2(v):
	return  np.mean(v)


def dstar_L1(v):
	return np.linalg.norm(v-np.median(v),1)

def mustar_L1(v):
	return  np.median(v)
