'''
Created on Feb 12, 2013

@author: lunt



I liberally copy from myself.
'''

from amaranth import simformat

from itertools import izip,count
import scipy as S
import scipy.sparse as sp


import Bio.Alphabet.IUPAC

alph = Bio.Alphabet.IUPAC.IUPACProtein.letters + "-"
Q = len(alph)
small_Q = Q-1
alphIndex = dict([(alph[i],i) for i in range(len(alph))])
alphIndex['X'] = alphIndex['-']
alphIndex['Z'] = alphIndex['-']
alphIndex['B'] = alphIndex['-']
def intConv(aString):
	return [alphIndex[i] for i in aString.strip()]

def weightedLoad(infile,weightthresh=None):
	
	myAlign = simformat.read(infile)
	myHeader = myAlign.header
	
	if weightthresh is not None:
		try:
			weightsindex = myHeader.cutoffs.index(weightthresh)
		except:
			raise Exception("No such weighting cutoff, valid cutoffs are: " + repr(myHeader.cutoffs))
		simformat.annotateAlignment(myAlign)
		weights = S.array(simformat.getinvnormsim(myAlign,weightsindex))
	else:
		weights = S.ones(len(myAlign))
	
	N = len(myAlign)
	Width = len(myAlign[0])

	Matrix = sp.lil_matrix((N,Q*Width)) #LiL is better to populate, csc might be even better, but we'd have to write more complex code.
	
	for seqRec,one_weight,i in izip(myAlign,weights,count()):
		seq_as_ints = intConv(seqRec.seq.tostring())
		for residue,j in izip(seq_as_ints,count()):
			Matrix[i,j*Q + residue] = one_weight

	return Matrix.tocsc()

def consensus(WeightedAlignMatrix,threshold=0.95,ambiguous='X'):
	PSSM = WeightedAlignMatrix.sum(0).A[0].reshape(-1,Q)
	TotalWeight = PSSM[0].sum()
	PSSM/=TotalWeight
	
	conservation = PSSM.max(1)
	alphabetIndeces = PSSM.argmax(1)
	consensus = ''.join([alph[i] if j >= threshold else ambiguous for i,j in izip(alphabetIndeces,conservation)])
	
	return consensus


if __name__ == "__main__":
	"""
	Do some basic tests
	"""
	
	import sys
	
	INFILE = open(sys.argv[1])
	
	AlignMat = weightedLoad(INFILE,0.7)
	
	for oneT in [0.95,0.9,0.85,0.8,0.7,0.5, 0.25, 0.2,0.1,0.05,0.03]:
		cons = consensus(AlignMat,threshold=oneT)
		print oneT, cons