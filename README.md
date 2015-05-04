EVFold Reimplementation

This is a clean-room reimplementation of EVFold ( www.evfold.org )


I am sorry about its state.


Running: (in the bash shell)
	1) Environment
		export PYTHONPATH=./utils/pymodules:${PYTHONPATH}
	2) Weighted (.sim) input file.
		./utils/SimpleHammWeight.py <input>.faa <output>.sim

		Be careful, on huge alignments this may take quite a lot of time and memory.
		Consider using hhblits to create your alignments and hhfilt to filter to get fewer sequences with more independence. That will speed things up a lot.
		My favorite settings for hhblits are " -n 4 -nodiff -cov 75 -id 99 -nodiff -neffmax 20 --maxfilt 100000 "
		My favorite settings for hhfilt are " hhfilter -cov 75 -id 95 -i <INTPUT> -o <OUTPUT>"
	3) 
		consider seeing bin/prepare.bash and Makefile to see how to run the program.
