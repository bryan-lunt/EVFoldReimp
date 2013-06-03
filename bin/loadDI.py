import scipy as S
import scipy.sparse as SP

def loadDI(INFILE):
	DI_matrix_by_index = S.loadtxt(INFILE).T
	size = DI_matrix_by_index.max(1)[1] + 1

	return SP.coo_matrix((DI_matrix_by_index[2],DI_matrix_by_index[:2]),shape=(size,size)).todense().A


if __name__ == "__main__":
	import sys
	import pylab
	
	from optparse import OptionParser
	parser = OptionParser()
	options,args = parser.parse_args()
	
	for filename in args:
		D = loadDI(filename)
		pylab.matshow(D)
		pylab.title(filename)

	pylab.show()