#!/bin/csh -f
  if($#argv != 4) then
echo "Usage: {nmashift} {input-pdb} {mode} {rms} {nbpas}" ; exit
  endif
set MODE = $2
set RMS = $3
set NMAX = $4
set BIN = $VEDA/$VEDA_BIN

cat >! DaTa <<EOF
matrix.eigenfacs
$1

n
y
$RMS PAS
$MODE
EOF

set n=-$NMAX
   while ($n <= $NMAX)
sed -e "s/PAS/${n}/" DaTa >! DaTo
ln -sf $PWD/out${n} fort.12
$BIN/projmod.exe < DaTo
set n=`expr $n + 1`
   end
ls -v out* >! DaTo

cat `cat DaTo` >! models.pdb

rm -f fort.12 DaTa DaTo

exit
