#!/bin/bash

var=$1
script=scripts/cntr-$var.gri
gpfile=$2
outfile=${gpfile%.gp}-$var.eps

timestamp=`echo $gpfile | tr -d 'a-z.' | sed "s/^0+//"|perl -p -e 's/^0+/0/;$_*=1e-5' \
	|perl -e 'while(<>) { printf("%.2f\n",$_);}'`
timestamp="t=$timestamp"
echo $timestamp $outfile
cat $script | sed "s/GPFILENAME/$gpfile/g;s/TIMESTAMP/$timestamp/g" > tmp.gri && \
	gri -output $outfile tmp.gri

