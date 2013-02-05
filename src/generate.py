#!/usr/bin/env

import scipy as S

import sys
from itertools import izip, count

import Bio.SeqIO
import Bio.PDB.Polypeptide as PP

from SS_const import *
from loadDI import *
from readSS import *


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

seqStr = seq.seq.tostring()
seqLen = len(seqStr)

outSeq = [PP.one_to_three(i) for i in seqStr]

OUTSEQ = open(OUTSEQ_FILENAME,'w')

for i in range(seqLen/10):
	OUTSEQ.write(' '.join(outSeq[i*10:i*10+10]))
	OUTSEQ.write('\n')
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
#TODO: Break into own function?
#TODO: see Table S2, this just masks all constraints within one SS element, but that is not exactly what they do in the paper.
for start, end, type in SS_map:
	if type is 'H':
		DI_Matrix[max(start-3,1):min(end+1,seqLen),max(start-3,1):min(end+1,seqLen)] = 0.0
	if type is 'E':
		DI_Matrix[max(start-4,1):min(end+2,seqLen),max(start-4,1):min(end+2,seqLen)] = 0.0


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

