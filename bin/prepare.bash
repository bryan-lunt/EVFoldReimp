#!/bin/bash

#####
#This script will do the basic preparation for a run.
#Predict the secondary structure,
#Create the enlongated sequence,
#Etc, Etc.
#
#Mostly, It wraps a python script.
#

SCRIPTDIR=`readlink -f $0`

CNSVERSION='1.3'

WORKDIR=.

####
#read command line arguments:
#TODO







####
#Set Environment variables that may depend on command line arguments.
#
#

CNSSCRIPTS=${SCRIPTDIR}/../cns_scripts/${CNSVERSION}/
SCRIPTOUTDIR=${WORKDIR}/output

####
#Confirm that expected files are present.

#single.faa -> sequence, secondary.ss, etc
#coupings.DI -> contacts
#alignment OR some_kind_of_conservation_file. -> masking of contacts


####
#Create some other directory structure
# ./output
#
mkdir -p ${SCRIPTOUTDIR}

####
#Predict SS
jnet -z single.faa > secondary.ss

####
#Generate constraints
${SCRIPTDIR}/generate.py single.faa couplings.DI secondary.ss alignment.sim
#
#Now, we must have the following files:
#	OUTSEQ_FILENAME = 'my.seq'
#	OUTSS_FILENAME = 'my_SS.tbl'
#	OUTSSA_FILENAME = 'my_SS_angle.tbl'
#	OUTCON_FILENAME = 'my.tbl'
#

####
#Generate an extended PDB file
cns_solve < ${CNSSCRIPTS}/generate_seq.inp > ${SCRIPTOUTDIR}/generate_seq.out
cns_solve < ${CNSSCRIPTS}/generate_extended.inp > ${SCRIPTOUTDIR}/generate_extended.out


####
#Test that expected files are present.
#TODO
#

