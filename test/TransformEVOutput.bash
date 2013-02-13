#!/bin/bash

INFILE=$1
sed ' /^!/ d; s/!.*//' ${INFILE} | sort -g -k8 | sort -g -k3 | tr '[a-z]' '[A-Z]'
