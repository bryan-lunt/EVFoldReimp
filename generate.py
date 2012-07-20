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
INFILE = open(sys.argv[1])

DI_Matrix = loadDI(sys.argv[2])

HMMTOPFILE = open(sys.argv[3])

#convert the FASTA sequence to input for CNS
seq = Bio.SeqIO.read(INFILE,'fasta')

seqStr = seq.seq.tostring()
seqLen = len(seqStr)

outSeq = [PP.one_to_three(i) for i in seqStr]

OUTSEQ = open('my.seq','w')

for i in range(seqLen/10):
	OUTSEQ.write(' '.join(outSeq[i*10:i*10+10]))
	OUTSEQ.write('\n')
OUTSEQ.close()
#end convert


#make SS constraints from HMMTOP
SSfile = open(sys.argv[3])
SS_map = make_ranges(read_jnet(SSfile))


#actually make SS constraints
OUTSS = open('my_SS.tbl','w')
OUTSSA = open('my_SS_angle.tbl','w')

for Hstart, Hend, type in SS_map:
	thisTM_dists = make_SS_dists(Hstart+1,Hend-Hstart+1)
	thisDihedral = make_helix_dihedral(Hstart+1,Hend-Hstart+1)

	OUTSS.write(thisTM_dists)
	OUTSSA.write(thisDihedral)

OUTSS.close()
OUTSSA.close()

#use SS map to mask DI constraints
for start, end, type in SS_map:
	if type in ['E','H']:
		DI_Matrix[start:end,start:end] = 0.0

nonZ = DI_Matrix.nonzero()
DIs = DI_Matrix[nonZ]

#make DI constraints
constraintstr = 'assign (resid %i and name CA)  (resid %i and name CA)  %.1f %.1f %.1f weight %.5f\n'

OUTCONS = open('my.tbl','w')
for i,j,k,number in izip(nonZ[0],nonZ[1],DIs,count(1)):
	OUTCONS.write(constraintstr % (i+1,j+1, 7.0, 4.0, 3.0,10.0/number))
OUTCONS.close()

