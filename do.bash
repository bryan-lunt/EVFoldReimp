#!/usr/bin/bash


cns_solve < cns_scripts/generate_seq.inp > cns_output/generate_seq.out
cns_solve < cns_scripts/generate_extended.inp > cns_output/generate_extended.out
cns_solve < cns_scripts/anneal.inp > cns_output/anneal.out
cns_solve < cns_scripts/accept.inp > cns_output/accept.out
cns_solve < cns_scripts/dg_sa.inp > cns_output/dg_sa.out
