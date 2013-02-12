prep:
	jnet -z foo.faa > foo.ss
	python src/generate.py foo.faa foo.DI foo.ss foo.sim
	cns_solve < cns_scripts/1.3/generate_seq.inp > cns_output/generate_seq.out
	cns_solve < cns_scripts/1.3/generate_extended.inp > cns_output/generate_extended.out

run:
	cns_solve < cns_scripts/1.3/anneal.inp > cns_output/anneal.out
	cns_solve < cns_scripts/1.3/accept.inp > cns_output/accept.out
	cns_solve < cns_scripts/1.3/dg_sa.inp > cns_output/dg_sa.out

all: clean prep run

clean:
	rm -f cns_output/* my* *.pdb *.pyc *.pyo
