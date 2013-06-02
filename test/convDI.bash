#!/bin/bash
#
#Convert a Sander/Marks MI_EC.txt file to our familiar DI file format. (Maybe?)
#
awk '{print $1, $3, $6}' $1
