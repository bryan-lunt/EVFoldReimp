#!/bin/bash
#
#Script to convert Marks/Sander et al XXXXX_ECScores.csv file into our familiar DCA file format.
#It is not yet clear if they are 0 or 1 indexed.
#
awk -F',' '{print $2, $1, $3}' $1 | sort -g -r -k3