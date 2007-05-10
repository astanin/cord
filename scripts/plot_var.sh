#!/bin/bash

if [ "foo$1" == "foo" ]; then
	script=scripts/gnuplot_phi
	echo "plotting phi"
else
	script=scripts/gnuplot_$1
fi

if [ -f log.dat ] ; then
	rm log.dat
fi

cord_hdf2gp *.hdf

for f in *.gp ; do 
	timestamp=`echo $f | tr -d 'a-z.' | sed "s/^0+//"`
	title="t=$timestamp"
	cat "$script" | sed "s/filename/$f/g;s/%title%/$title/" | gnuplot - 2>>log.dat
done

