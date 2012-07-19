#!/usr/bin/env

import scipy as S

import sys
from itertools import izip, count

import Bio.SeqIO
import Bio.PDB.Polypeptide as PP

from SS_const import *


INFILE = open(sys.argv[1])

DI = S.loadtxt(sys.argv[2])
DI[:,0:2] += 1

seq = Bio.SeqIO.read(INFILE,'fasta')

seqStr = seq.seq.tostring()
seqLen = len(seqStr)


outSeq = [PP.one_to_three(i) for i in seqStr]

OUTSEQ = open('my.seq','w')

for i in range(seqLen/10):
	OUTSEQ.write(' '.join(outSeq[i*10:i*10+10]))
	OUTSEQ.write('\n')
OUTSEQ.close()


#make DI constraints
constraintstr = 'assign (resid %i and name CA)  (resid %i and name CA)  %.1f %.1f %.1f weight %.5f\n'


OUTCONS = open('my.tbl','w')
for (i,j,k),number in izip(DI,count(1)):
	OUTCONS.write(constraintstr % (int(i),int(j), 7.0, 4.0, 3.0,10.0/number))
OUTCONS.close()


#make SS constraints from HMMTOP
HMMTOPFILE = open(sys.argv[3])
for line in HMMTOPFILE:
	if line.startswith('Transmembrane helices'):
		TMstr = map(lambda x: tuple(map(int,x.split('-'))), line.split(':')[1].strip().split(' '))
HMMTOPFILE.close()


#actually make SS constraints
OUTSS = open('my_SS.tbl','w')
OUTSSA = open('my_SS_angle.tbl','w')

for Hstart, Hend in TMstr:
	thisTM_dists = make_SS_dists(Hstart,Hend-Hstart)
	thisDihedral = make_helix_dihedral(Hstart,Hend-Hstart)

	OUTSS.write(thisTM_dists)
	OUTSSA.write(thisDihedral)

OUTSS.close()
OUTSSA.close()
