#!/usr/bin/env

import scipy as S

import sys
from itertools import izip, count, groupby

import consensus

import Bio.SeqIO
import Bio.PDB.Polypeptide as PP

from contact_const import *
from SS_const import *
from loadDI import *
from readSS import *


def fastaToCNSThreeLetter(seqRecord,width=10):
	seqStr = seqRecord.seq.tostring()
	L = len(seqStr)
	
	outSeq = [PP.one_to_three(i) for i in seqStr]
	numbered = zip(outSeq,count(0))
	
	groups = [[j[0] for j in i[1]] for i in groupby(numbered,lambda x:x[1]/width)]
	joined = '\n'.join([' '.join(one_group) for one_group in groups]) + '\n'
	
	return joined
	


def generateDIMask(theDIs, cons, sStruct, k=6,consthresh=0.95,ambigChar='X'):
	L,M = theDIs.shape
	assert L == M , "Non square DI matrix!?"
	DImask = S.ones(theDIs.shape)
	
	#
	#1 Sequence Separation
	#
	DImask = S.triu(DImask,k)
	
	#
	#2 Conservation Filter
	#
	#Scan the consensus sequence and find overly conserved columns, note those, and note conserved cysteine columns
	conservedSites = S.array([i != ambigChar for i in cons],dtype=S.int8)
	conservedCysteine = S.array([i == 'C' for i in cons],dtype=S.int8)
	
	DImask = (1-conservedSites).reshape(-1,1)*(DImask*(1-conservedSites))
	
	#
	#3 Cysteine columns may be strongly conserved, but only allow their best contact.
	#TODO : implement Cysteine Best-Friend pairs
	
	
	#
	#4 Secondary Structure Masking
	#Seet Table S2 in th EVFold paper
	#
	abovez = lambda x:max(x,0)
	for start,end,type in sStruct:
		if type is "H":
			#alpha helix
			#Table S2 in EVFold paper. Not how I would have done it.
			DImask[abovez(start-1):end+1,abovez(start-1):end+1] = 0.0 #Constraints H,H; H-1,H; H-1,H+1 in Table S2
			DImask[abovez(start-2):abovez(end-4),start+4:end+2] = 0.0 #Constraints H-2,H; (assumed) H,H+2 in Table S2
			DImask[abovez(start-3):abovez(end-9),start+9:end+3] = 0.0 #Constraints H-3,H; (assumed) H,H+3 in Table S2
			#TODO : Not mentioned in the paper, but what if helices are longer than this? We could just continue in a pattern!
		elif type is "E":
			#beta strand
			width = 1 if (end-start) < 5 else 2
			DImask[abovez(start-width):end+width,abovez(start-width):end+width] = 0.0 #constraints E,E ; E-1, E+1; E-2,E+2 in Table S2
			DImask[abovez(start-2):end,start:end] = 0.0#constraint E-2,E
			DImask[start:end,start:end+2] = 0.0#constraints (assumed) E,E+2
			DImask[abovez(start-3):end,start:end] = 0.0#constraint E-3,E
			DImask[start:end,start:end+3] = 0.0#constraints (assumed) E,E+3
			DImask[abovez(start-4):start,start+4:end] = 0.0#constraint E-4,E in table S2
			DImask[start:abovez(end-4),end:end+4] = 0.0#constraint (assumed) E,E+4 in table S2i
	
	
	return DImask
	

from optparse import OptionParser
def main():
	parser = OptionParser(usage='usage: %prog [options] <SingleSequence.faa> <DI file> <SS Prediction> <SIMformat Alignment>')
	
	parser.add_option("--wt","--weight-thresh",dest="weightthresh",help="The weighting threshold to use if using reweighting.",default=0.7,type="float",metavar="FLOAT")
	parser.add_option("--ct","--consensus-thresh",dest="consthresh",help="The conservation threshold at which DI pairs for that column are ignored (unless cysteine)",default=0.95,type="float",metavar="FLOAT")
	parser.add_option("-n","--number",dest='numpairs',help='The number of pairs to take',type='int',default=None)
	parser.add_option('-f','--fract',dest='fraction',help='The number of pairs will be Fract*L where Fract is this option and L is the length of the sequence',type='float',default=0.2)
	
	options, args = parser.parse_args()
	
	#list input files
	INSEQ = open(args[0])
	DI_Matrix = loadDI(args[1])
	M,L = DI_Matrix.shape
	SSPREDFILE = open(args[2])
	INALIGN = open(args[3])
	
	WeightThresh = options.weightthresh
	ConsensusThresh = options.consthresh
	NumPairs = options.numpairs if options.numpairs is not None else int(round(options.fraction*L)) #TODO: Calculate this from L, the length of the sequence
	
	OUTSEQ_FILENAME = 'my.seq'
	OUTSS_FILENAME = 'my_SS.tbl'
	OUTSSA_FILENAME = 'my_SS_angle.tbl'
	OUTCON_FILENAME = 'my.tbl'
	#
	#General Input
	#
	seq = Bio.SeqIO.read(INSEQ,'fasta')
	align = consensus.weightedLoad(INALIGN,WeightThresh)
	cons = consensus.consensus(align,ConsensusThresh)
	#read SS from HMMTOP/JNET
	#TODO: Commandline option to use HMMTOP/JNET or autodetect...
	SS_seq = read_jnet(SSPREDFILE)
	SS_map = make_ranges(SS_seq)
	
	
	
	
	#convert the FASTA sequence to input for CNS
	outSequence = fastaToCNSThreeLetter(seq)
	OUTSEQ = open(OUTSEQ_FILENAME,'w')
	OUTSEQ.write(outSequence)
	OUTSEQ.close()
	#end convert

	
	
	
	#actually make SS constraints
	OUTSS = open(OUTSS_FILENAME,'w')
	OUTSSA = open(OUTSSA_FILENAME,'w')
	
	for Hstart, Hend, type in SS_map:
		thisTM_dists = make_SS_dists(Hstart+1,Hend-Hstart+1,type=type)
		OUTSS.write(thisTM_dists)
		
		thisDihedral = make_SS_angles(Hstart+1,Hend-Hstart+1,type=type)
		OUTSSA.write(thisDihedral)
	
	OUTSS.close()
	OUTSSA.close()
	
	
	
	#use SS map to mask DI constraints
	mask = generateDIMask(DI_Matrix,cons,SS_map)
	print "%i positions survived masking." % mask.sum()
	DI_Matrix = S.multiply(DI_Matrix,mask)#Hadamard Product!
	
	M,L = DI_Matrix.shape
	assert M == L, 'DI matrix not square!?'
	DI_Ravel = DI_Matrix.ravel()
	
	print "Taking %i pairs" % NumPairs
	
	topIndeces = DI_Ravel.argsort()[-1:0:-1][:NumPairs]
	topPositions = S.array([(i/L,i%L) for i in topIndeces],dtype=int)
	
	pair_constraints = make_EVFold_contacts(topPositions,seq.seq.tostring())
	
	OUTCONS = open(OUTCON_FILENAME,'w')
	OUTCONS.write(pair_constraints)
	OUTCONS.close()

if __name__ == "__main__":
	main()
