#!/bin/bash

INFILE=$1
# sed deletes any comments or commented lines in the output file
# sort puts them in the same order
# tr transforms them to upper case
sed ' /^!/ d; s/!.*//' ${INFILE} | sort -g -k8 | sort -g -k3 | tr '[a-z]' '[A-Z]'
