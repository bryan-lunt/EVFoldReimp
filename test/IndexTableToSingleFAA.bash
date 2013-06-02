#!/bin/bash
#
#Generate a single-sequence .faa file from the Marks/Sander XXXX.indextable file.
#
sed '1 d' $1 | awk 'START { foo = ""; }; {foo = foo  $2}; END {print ">seq"; print foo;}'
