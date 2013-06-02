#!/bin/bash
#
#Recover Jnet SS file from the Marks/Sander XXXX.indextable file.
#
sed '1 d' $1 | awk 'START { print "START PRED";}; {bar = $3; if (bar == "C") bar = "-"; print $2 " " bar " | " bar " " ;}; END {print "END PRED"}'
