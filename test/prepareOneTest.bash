#!/bin/bash
#
#This should prepare and execute one regression test.
#

#The directory containing a bundle of EVFold supplementary output.
INDIR=$1
OUTDIR=$2

#preparing input
#Alignment: Strip the a2m file, call the hamming reweighter.
stripInserts.sed ${INDIR}/*a2m > ${OUTDIR}/stripped.faa
SimpleHammWeight.py  ${OUTDIR}/stripped.faa  ${OUTDIR}/stripped.sim

#prepare the DI
awk '{print $1, $3, $6}' ${INDIR}/*MI_ECs.txt | sort -g -r -k3 > ${OUTDIR}/theDI.DI

#prepare the secondary structure
sed '1 d' ${INDIR}/*.indextable | awk 'START { print "START PRED";}; {bar = $3; if (bar == "C") bar = "-"; print $2 " " bar " | " bar " " ;}; END {print "END PRED"}' > ${OUTDIR}/theSS.SS

#prepare the single file
sed '1 d' ${INDIR}/*.indextable | awk 'START { foo = ""; }; {foo = foo  $2}; END {print ">seq"; print foo;}' >  ${OUTDIR}/theFAA.faa