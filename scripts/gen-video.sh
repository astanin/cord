#!/bin/bash

if [ $# -le 1 ]; then
	echo usage:
	echo gen-video.sh frames-per-second frame1.img \[frame2.img ... \]
	exit -1
fi

fps=$1
frame=0;
shift

while [ $# -gt 0 ]; do
	timestamp=""
	printf -v timestamp "%08d" $frame
	framename="frame.$timestamp.jpg"
	convert -quality 95 "$1" "$framename"
	echo $framename: $1
	shift
	frame=$((frame+1));
done

mencoder "mf://frame.*.jpg" -mf fps=$fps -zoom -xy 1 -oac copy -ovc lavc \
	-lavcopts vcodec=msmpeg4:vbitrate=4000 -ffourcc MP43 -forceidx \
	-o output.avi

rm -rf frame.[0-9]*.jpg
