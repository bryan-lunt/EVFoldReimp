#!/usr/bin/env python
"""
Use this program to calculate the weight of faa alignements based on all vs all hamming distances.

Usage: %prog [options] <inputFAA> <outputSIM>
"""
import sys
from itertools import izip, imap, count
from optparse import OptionParser

import scipy as S
import scipy.sparse as sp
import numpy as N

import scipy.spatial.distance as D

import amaranth.simformat as simformat

def uniqSeqs(aList,key=lambda x: x.seq.tostring()):
	"""
	This function will return a list of only the unique sequence records, by whatever key the user supplies, usually the sequence itself.
	"""
	seenSeqs = set()
	def haveSeen(aSeq):
		if aSeq in seenSeqs:
			return True
		else:
			seenSeqs.add(aSeq)
			return False
	foo = [i for i in aList if not haveSeen(key(i))]
	return foo



#main
def main():
	#parse commandline arguments
	parser = OptionParser(usage="Usage: %prog [options] <inputFAA> <outputSIM>")
	parser.add_option('-u','--unique',dest='unique',default=True,action='store_true')
	parser.add_option('-i','--ids',dest='cutoffs',type='string',default='100,98,95,90,85,80,75,70')
	parser.add_option('-f','--fids',dest='fcutoffs',type='string',default=None)
	options, args = parser.parse_args()
	if len(args) != 2:
		parser.print_help()
		sys.exit()
		
	#read the alignment from file
	with open(args[0]) as infile:
		myAlignment = simformat.read(infile)

	if options.unique:
		myAlignment._records = uniqSeqs(myAlignment._records)
	
	simformat.annotateAlignment(myAlignment)
	
	if options.fcutoffs is None:
		thresholds = S.array(map(float,options.cutoffs.split(',')))/100.0
	else:
		thresholds = S.array(map(float,options.fcutoffs.split(',')))
	myAlignment.header.cutoffs = thresholds
	
	#MORE CODE HERE!!!
	AsVects = [S.array(map(ord, record.seq.tostring())) for record in myAlignment]
	
	AllSimilarities = 1.0 - D.cdist(AsVects,AsVects,'hamming')
	
	weights = S.zeros((len(myAlignment),len(thresholds)))
	
	for oneThresh, col in izip(thresholds,count()):
		weights[:,col] = (AllSimilarities > oneThresh).sum(1)
	
	for i in range(0,len(myAlignment)):
		myAlignment._records[i].annotations["weights"] = map(int,weights[i])
	
	#write output file
	with open(args[1],"w") as outfile:
		simformat.write(outfile,myAlignment)

if __name__ == "__main__":
	main()