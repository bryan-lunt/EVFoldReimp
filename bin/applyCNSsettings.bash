#!/bin/bash
#

#Use this script to replace variables in CNS scripts
#
#You should set the varible "FOLD_SUBJOB"
#

BLUNT_GEN_SEQ_INFILE="my.seq"
BLUNT_GEN_SEQ_OUTFILE="my.mtf"

BLUNT_GEN_EXT_OUTFILE="my_extended.pbd"

BLUNT_CONSTRAINTS_TABLE="my.tbl"
BLUNT_SS_CONSTRAINTS_TABLE="my_SS.tbl"
BLUNT_SS_ANGLES_TABLE="my_SS_angle.tbl"

#for DG_SA Distance Geometry Embedding + Simulated Annealing
BLUNT_EMBED_HOWMANY_PER_TASK="5"
BLUNT_EMBEDDED_ANNEALED_BASE="${FOLD_SUBJOB}_dg"

#for Straight simulated annealing
BLUNT_ANNEAL_BASE="${FOLD_SUBJOB}_anneal"
BLUNT_ACCEPT_BASE="${FOLD_SUBJOB}_accept"

#for minimization in original (EVFOLD) files
BLUNT_MIN_BASE="${FOLD_SUBJOB}_min"


###
#Make Substitutions
#
sed "s/BLUNT_GEN_SEQ_INFILE/${BLUNT_GEN_SEQ_INFILE}/g; \
s/BLUNT_GEN_SEQ_OUTFILE/${BLUNT_GEN_SEQ_OUTFILE}/g; \
s/BLUNT_GEN_EXT_OUTFILE/${BLUNT_GEN_EXT_OUTFILE}/g; \
s/BLUNT_CONSTRAINTS_TABLE/${BLUNT_CONSTRAINTS_TABLE}/g; \
s/BLUNT_SS_CONSTRAINTS_TABLE/${BLUNT_SS_CONSTRAINTS_TABLE}/g; \
s/BLUNT_SS_ANGLES_TABLE/${BLUNT_SS_ANGLES_TABLE}/g; \
s/BLUNT_EMBED_HOWMANY_PER_TASK/${BLUNT_EMBED_HOWMANY_PER_TASK}/g; \
s/BLUNT_EMBEDDED_ANNEALED_BASE/${BLUNT_EMBEDDED_ANNEALED_BASE}/g; \
s/BLUNT_ANNEAL_BASE/${BLUNT_ANNEAL_BASE}/g; \
s/BLUNT_ACCEPT_BASE/${BLUNT_ACCEPT_BASE}/g; \
" $1