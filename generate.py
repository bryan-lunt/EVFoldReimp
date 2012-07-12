#!/usr/bin/env

import scipy as S

import sys

import Bio.SeqIO
import Bio.PDB.Polypeptide as PP

INFILE = open(sys.argv[1])

DI = S.loadtxt(sys.argv[2])
DI[:,0:2] += 1

seq = Bio.SeqIO.read(INFILE,'fasta')

outSeq = ' '.join([PP.one_to_three(i) for i in seq.seq.tostring()])

OUTSEQ = open('my.seq','w')
OUTSEQ.write(outSeq)
OUTSEQ.write('\n')
OUTSEQ.close()


#make DI constraints
constraintstr = 'assign (resid %i and name CA)  (resid %i and name CA)  %.1f %.1f %.1f \n'


OUTCONS = open('my.tbl','w')
for i,j,k in DI:
	OUTCONS.write(constraintstr % (int(i),int(j), 7.0, 1.0,2.0))
OUTCONS.close()


#make SS constraints from HMMTOP
HMMTOPFILE = open(sys.argv[3])
for line in HMMTOPFILE:
	if line.startswith('Transmembrane helices'):
		TMstr = map(lambda x: tuple(map(int,x.split('-'))), line.split(':')[1].strip().split(' '))
HMMTOPFILE.close()

OUTSS = open('my_SS.tbl','w')
for Hstart, Hend in TMstr:
	for i in range(Hstart,Hend+1):
		if i+1 <= Hend:
			OUTSS.write(constraintstr % (i, i+1, 3.82, 0.2, 0.2))
		if i+2 <= Hend:
			OUTSS.write(constraintstr % (i, i+2, 5.5, 0.3, 0.3))
		if i+3 <= Hend:
			OUTSS.write(constraintstr % (i, i+3, 5.3, 0.65, 0.65))
		if i+4 <= Hend:
			OUTSS.write(constraintstr % (i, i+4, 6.35, 0.7, 0.7))
		if i+5 <= Hend:
			OUTSS.write(constraintstr % (i, i+5, 8.7, 0.6, 0.6))
OUTSS.close()
