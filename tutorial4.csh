#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/trx/trx.mtf \
	     data/trx/trx_extended.pdb \
	     data/trx/trx_noe_all.tbl \
	     data/trx/trx_noe_hbond.tbl \
	     data/trx/trx_dihed_rama.tbl \
             nmr/dgsa/dg_sa.inp)
#
if ( ! $?CNS_SOLVE ) then
  echo "CNS_SOLVE not defined"
  echo "CNS must be correctly setup to running this tutorial script"
  exit 1
endif
#
if ( ! -d $CNS_SOLVE/$tutdir ) then
  echo "cannot find tutorial directory $CNS_SOLVE/$tutdir"
  echo "please check that you are using the correct version of CNS"
  exit 1
endif
#
foreach file ( $lnfiles )
  if ( ! -e {$file}:t ) then
    echo "linking $tutdir/$file to current directory"
    ln -s $CNS_SOLVE/$tutdir/$file .
  endif
end
#
echo "        "
echo "running:"
echo "         cns_solve < dg_sa.inp > dg_sa.out"
echo "        "
#
cns_solve < dg_sa.inp > dg_sa.out
#
echo "        "
#
