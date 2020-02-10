# EVFold Reimplementation

This is an unfinished clean-room reimplementation of EVFold ( www.evfold.org )

Unfortunately Sanders/Marks distribute the DCA code [Martin Weigt's scientific contribution], but not the folding code that works on top of DCA output [their scientific contribution]. The only way to test if their method even works as described has been to try to reverse-engineer/reimplement it from what few SI files were made available. Recently they have taken down those files, but you can still find them on www.archive.org .

I am sorry about the state of this. If people really worked in an attitude of openness, it would never have been necessary in the firstplace.

## Installation

	You need Python to run the scripts. They are old and may only work on py2.7 . (sorry again)

	You need the CNS package to do the computation. Licensing prevents me from providing it to you. [ https://www.mrc-lmb.cam.ac.uk/public/xtal/doc/cns/cns_1.3/installation/frame.html ]
	There are scripts for CNS version 1.3 and 1.21 . I hope they both work.
	Again, I would have liked to provide a docker image for this, but licensing prevents my redistributing CNS. :/



## Running: (in the bash shell)

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
