#!/bin/bash

sed '1 d' $1 | awk 'START { foo = ""; }; {foo = foo  $2}; END {print ">seq"; print foo;}'
