#!/usr/bin/env

import scipy as S

import sys
from itertools import izip, count, groupby

import Bio.SeqIO
import Bio.AlignIO
from Bio.Align import AlignInfo
import Bio.PDB.Polypeptide as PP

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
	


def generateDIMask(theDIs, theAlign, sStruct, k=6,consthresh=0.95):
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
	myInfo = AlignInfo.SummaryInfo(theAlign)
	cons = myInfo.gap_consensus(threshold = consthresh,ambiguous='X')
	#Scan the consensus sequence and find overly conserved columns, note those, and note conserved cystein columns
	conservedSites = S.zeros(L)
	conservedCysteine = S.zeros(L)
	for i in range(L):
		res = cons.seq[i]
		if res is not "X":
			if res is "C":
				conservedCysteine[i] = 1.0
			else:
				conservedSites[i] = 1.0
	
	allConservation = conservedSites + conservedCysteine
	
	DImask = allConservation.reshape(-1,1)*(DImask*allConservation)
	
	
	#
	#3 Cysteine columns may be strongly conserved, but only allow their best contact.
	#TODO : implement this later.
	
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
	

def main():
	#list input files
	INSEQ = open(sys.argv[1])
	DI_Matrix = loadDI(sys.argv[2])
	SSPREDFILE = open(sys.argv[3])
	
	OUTSEQ_FILENAME = 'my.seq'
	OUTSS_FILENAME = 'my_SS.tbl'
	OUTSSA_FILENAME = 'my_SS_angle.tbl'
	OUTCON_FILENAME = 'my.tbl'
	
	#convert the FASTA sequence to input for CNS
	seq = Bio.SeqIO.read(INSEQ,'fasta')
	
	outSequence = fastaToCNSThreeLetter(seq)
	
	OUTSEQ = open(OUTSEQ_FILENAME,'w')
	OUTSEQ.write(outSequence)
	OUTSEQ.close()
	#end convert

	
	
	#make SS constraints from HMMTOP/JNET
	#TODO: Commandline option to use HMMTOP/JNET or autodetect...
	SS_map = make_ranges(read_jnet(SSPREDFILE))
	
	
	#actually make SS constraints
	OUTSS = open(OUTSS_FILENAME,'w')
	OUTSSA = open(OUTSSA_FILENAME,'w')
	
	for Hstart, Hend, type in SS_map:
		thisTM_dists = make_SS_dists(Hstart+1,Hend-Hstart+1,type=type)
		thisDihedral = make_SS_angles(Hstart+1,Hend-Hstart+1,type=type)
	
		OUTSS.write(thisTM_dists)
		OUTSSA.write(thisDihedral)
	
	OUTSS.close()
	OUTSSA.close()
	
	#use SS map to mask DI constraints
	mask = generateDIMask(DI_Matrix,None,SS_map)

	
	
	#make DI constraints
	#TODO: Break this out into its own function, as it is starting to get complex.
	#TODO: See the paper supplementary text page 9, we need a conservation filter, that allows disulphide bonds despite high conservation, but only allows zero or one disulphide bonds per column
	
	
	#TODO: We also need a constraint between CB in pairs that do not contain glycine.
	constraintstr = 'assign (resid %i and name CA)  (resid %i and name CA)  %.1f %.1f %.1f weight %.5f\n'
	nonZ = DI_Matrix.nonzero()
	
	OUTCONS = open(OUTCON_FILENAME,'w')
	for i,j,number in izip(nonZ[0],nonZ[1],count(1)):
		k = DI_Matrix[i,j]
		OUTCONS.write(constraintstr % (i+1,j+1, 4.0, 4.0, 3.0, 10.0/number))
	OUTCONS.close()

if __name__ == "__main__":
	main()
