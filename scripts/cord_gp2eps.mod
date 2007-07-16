#!/bin/sh

PLOTTER=gri
MODE=mono
TITLE=no
export PLOTHEIGHT=10
export PLOTWIDTH=10

function usage {
	echo usage:
	echo     $(echo $0|sed 's/.*\///g') \[-gri\|-gnuplot\] \[-mono\|-color\] \[-withtitle\] \[-xrange xmin:xmax\] \[-yrange ymin:ymax\] datafile1.gp \[ datafile2.gp ... \]
	exit -1
}

if [ $# -lt 1 ]; then
	usage
fi

function plot_c {
	f=$1;
	fo=$(echo $TIME|perl -e 'while(<>){printf("Oxygen_%07.0f.eps",10*$_);}');
	gri -b -c1 -output "$fo" << ENDGRI
#set x size $(echo "scale=2;$PLOTHEIGHT*($XMAX-$XMIN)/($YMAX-$YMIN)"|bc -s)
#set y size $PLOTHEIGHT
set x size $PLOTWIDTH
set y size $(echo "scale=2;$PLOTWIDTH*($YMAX-$YMIN)/($XMAX-$XMIN)"|bc -s)
set axes style 1
set axes style offset 0.5
set x axis $XMIN $XMAX 0.1 0.02 # $(perl -e "printf('%g',($XMAX-$XMIN)*0.1);")
set y axis $YMIN $YMAX 0.1 0.02 # $(perl -e "printf('%g',($YMAX-$YMIN)*0.1);")
draw axes
# c
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NC],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
${PLOT_C_DATA}
draw axes frame
# psi
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NPSI],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
set line width rapidograph 1
draw contour 0.0 unlabelled
${DRAW_C_TITLE}
ENDGRI
}

function plot_phi {
	f=$1
	fo=$(echo $TIME|perl -e 'while(<>){printf("Phi_%07.0f.eps",10*$_);}');
	gri -b -c0 -output "$fo" << ENDGRI
#set x size $(echo "scale=2;$PLOTHEIGHT*($XMAX-$XMIN)/($YMAX-$YMIN)"|bc -s)
#set y size $PLOTHEIGHT
set x size $PLOTWIDTH
set y size $(echo "scale=2;$PLOTWIDTH*($YMAX-$YMIN)/($XMAX-$XMIN)"|bc -s)
set axes style 1
set axes style offset 0.5
set x axis $XMIN $XMAX 0.1 0.02 # $(perl -e "printf('%g',($XMAX-$XMIN)*0.1);")
set y axis $YMIN $YMAX 0.1 0.02 # $(perl -e "printf('%g',($YMAX-$YMIN)*0.1);")
draw axes
# phi
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NPHI],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
${PLOT_PHI_DATA}
draw axes frame
# psi
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NPSI],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
set line width rapidograph 1
draw contour 0.0 unlabelled
${DRAW_PHI_TITLE}
ENDGRI
}

function plot_phi_growth {
	f=$1
	fo=$(echo $TIME|perl -e 'while(<>){printf("Growth_%07.0f.eps",10*$_);}');
	gri -b -c0 -output "$fo" << ENDGRI
#set x size $(echo "scale=2;$PLOTHEIGHT*($XMAX-$XMIN)/($YMAX-$YMIN)"|bc -s)
#set y size $PLOTHEIGHT
set x size $PLOTWIDTH
set y size $(echo "scale=2;$PLOTWIDTH*($YMAX-$YMIN)/($XMAX-$XMIN)"|bc -s)
set x axis $XMIN $XMAX $(perl -e "printf('%g',($XMAX-$XMIN)*0.1);")
set y axis $YMIN $YMAX $(perl -e "printf('%g',($YMAX-$YMIN)*0.1);")
draw axes
# phi
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NGROWTH],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
${PLOT_GROWTH_DATA}
# psi
set x grid $XMIN $XMAX /$XSIZE
set y grid $YMIN $YMAX /$YSIZE
open "cat $f|grep -v ^#|perl -e 'while(<>){if (\$_ !~ m/^$/ ) { @v=split; print \$v[$NPSI],\" \"; } else { print \"\n\";} }'|"
read grid data $XSIZE $YSIZE bycolumns
close
set line width rapidograph 1
draw contour 0.0 unlabelled
${DRAW_GROWTH_TITLE}
ENDGRI
}


function gri_plot {
	f=$1
	if [ "x$MODE" == "xmono" ]; then
	PLOT_C_DATA=$(cat << END
set image range 0.3 1.0
convert grid to image
set image greyscale black 0.3 white 0.9 increment 0.1
set image missing value color to rgb 1 0 0
draw image
draw image palette axisright left 0.4 right 1.0 increment 0.1 box 16.5 6.0 17.0 16.0 #box 3.5 6.0 4.0 16.0
set font size 9
set line width rapidograph 6x0
draw contour 0.4 1.0 0.1
END);
	PLOT_PHI_DATA=$(cat << END
set image range 0.67 0.83
convert grid to image
#set image greyscale black 0.69 white 0.81 increment 0.015
set image greyscale black 0.83 white 0.67 increment 0.02
set image missing value color to black
draw image
#draw image palette axisright left 0.705 right 0.795 increment 0.015 box 16.5 6.0 17.0 16.0 #box 3.5 6.0 4.0 16.0
draw image palette axisright left 0.81 right 0.69 increment 0.02 box 16.5 6.0 17.0 16.0 #box 3.5 6.0 4.0 16.0
set font size 9
set line width rapidograph 6x0
#draw contour 0.705 0.795 0.015
draw contour 0.69 0.81 0.02
END);
	PLOT_GROWTH_DATA=$(cat << END
#set font size 9
#set line width rapidograph 6x0
#draw contour
set image range -0.01 0.01
convert grid to image
open "LC_ALL=C awk 'BEGIN { for(i=0;i<256;i++) { if (i < (128-20)/2) { print 0.0; } else if (i < (128-20)) { print 0.15; } else if (((256-i) > (128-20)/2) && (i > (128+20))) { print 0.85;} else if (i > (128+20)) { print 0.7; } else { print 1.0; } } }' |"
read image greyscale
close
draw image
draw image palette axisleft left -0.01 right 0.01 box 3.5 6 4.5 16
END);
	else # colormode
	PLOT_C_DATA=$(cat << END
set image range 0.0 1.0
convert grid to image
set image colorscale rgb 0.2 0.2 0.7 1.00 rgb 1.0 1.0 0.5 0.0 increment 0.025
set image missing value color to rgb 1 0 0
draw image
draw image palette axisleft left 0.0 right 1.0 box 3.5 6.0 4.5 16.0
END);
	PLOT_PHI_DATA=$(cat << END
set image range 0.6 0.9
convert grid to image
set image colorscale rgb 1.0 0.894 0.882 0.6 rgb 0.604 0.161 0.161 0.9 increment 0.025
set image missing value color to rgb 0 1 0
draw image
draw image palette axisleft left 0.6 right 0.9 box 3.5 6 4.5 16
END);
	PLOT_GROWTH_DATA=$(cat << END
set image range -0.06 0.06
convert grid to image
#set image colorscale rgb 0.6 0.0 0.0 -0.06 rgb 0.7 0.7 1.0 0.06 increment 0.01
open "LC_ALL=C awk 'BEGIN { for(i=0;i<256;i++) { print(1-i/(3*256),(i-128)*(i-128)/(128*128),1.0+0.0*i/256) } }' |"
read image colorscale hsb
draw image
draw image palette axisleft left -0.06 right 0.06 box 3.5 6 4.5 16
END);
	fi
	if [ "x$TITLE" != "xno" ]; then
		TS=$(echo $TIME|perl -e "while(<>){printf('t=%5.1f',\$_);}");
		DRAW_C_TITLE=$(cat << END
set font size 12
draw title "Oxygen, $TS"
END);
		DRAW_PHI_TITLE=$(cat << END
set font size 12
draw title "Packing density, $TS"
END);
		DRAW_GROWTH_TITLE=$(cat << END
set font size 12
draw title "Rate of growth, $TS"
END);
	else
		DRAW_C_TITLE=""
		DRAW_PHI_TITLE=""
		DRAW_GROWTH_TITLE=""
	fi
	plot_c $f;
	plot_phi $f;
	plot_phi_growth $f;
}

function gnuplot_plot_c {
	f=$1
	fo=$(echo $TIME|perl -e 'while(<>){printf("Oxygen_%07.0f.eps",10*$_);}');
	gnuplot << END
set view 45,135
set xlabel "x"
set ylabel "y"
unset zlabel
unset key
${PLOT_C_TITLE}
set output "$fo"
${PLOT_MODE}
splot "$f" u 1:2:$((NC+1)) w l
END
}

function gnuplot_plot_phi {
	f=$1
	fo=$(echo $TIME|perl -e 'while(<>){printf("Phi_%07.0f.eps",10*$_);}');
	gnuplot << END
set view 45,135
set xlabel "x"
set ylabel "y"
unset zlabel
unset key
${PLOT_PHI_TITLE}
set output "$fo"
${PLOT_MODE}
splot "$f" u 1:2:$((NPHI+1)) w l
END
}

function gnuplot_plot_psi {
	f=$1
	fo=$(echo $TIME|perl -e 'while(<>){printf("Phase_%07.0f.eps",10*$_);}');
	gnuplot << END
set view 45,135
set xlabel "x"
set ylabel "y"
unset zlabel
unset key
${PLOT_PSI_TITLE}
set output "$fo"
${PLOT_MODE}
splot "$f" u 1:2:$((NPSI+1)) w l
END
}

function gnuplot_plot {
	f=$1
	if [ "x$MODE" == "xcolor" ]; then
		PLOT_MODE="set terminal post eps 14 color"
	else
		PLOT_MODE="set terminal post eps 14"
	fi
	if [ "x$TITLE" != "xno" ]; then
		TS=$(echo $TIME|perl -e "while(<>){printf('t=%5.1f',\$_);}");
		PLOT_C_TITLE="set title \"Oxygen, $TS\""
		PLOT_PHI_TITLE="set title \"Packing density, $TS\""
		PLOT_PSI_TITLE="set title \"Level set function, $TS\""
	else
		PLOT_C_TITLE=""
		PLOT_PHI_TITLE=""
		PLOT_PSI_TITLE=""
	fi
	gnuplot_plot_c $f
	gnuplot_plot_phi $f
	gnuplot_plot_psi $f
}

while [ $# -gt 0 ]; do
	f=$1;
	shift;
	opt=$f;
	if [ "x${opt:0:1}" == "x-" ]; then
		if [ "x$opt" == "x-color" ]; then
			MODE=color
			continue;
		elif [ "x$opt" == "x-mono" ]; then
			MODE=mono
			continue;
		elif [ "x$opt" == "x-withtitle" ]; then
			TITLE=yes
			continue;
		elif [ "x$opt" == "x-gnuplot" ]; then
			PLOTTER=gnuplot
			continue;
		elif [ "x$opt" == "x-gri" ]; then
			PLOTTER=gri
			continue;
		elif [ "x$opt" == "x-xrange" ]; then
			range=$1
			shift
			XMIN=${range%:*}
			XMAX=${range/*:/}
			CUSTOM_XRANGE=Y
			continue;
		elif [ "x$opt" == "x-yrange" ]; then
			range=$1
			shift
			YMIN=${range%:*}
			YMAX=${range/*:/}
			CUSTOM_YRANGE=Y
			continue;
		else
			echo "unknown option: $opt"
			usage
			break
		fi
	fi
	TIME=$(head -1 $f|tail -1|LC_ALL=C awk '{print $2}'|\
		perl -e 'while(<>){printf("%g\n",$_);}');
	XSIZE=$(head -2 $f|tail -1|LC_ALL=C awk '{print $2}');
	YSIZE=$(head -2 $f|tail -1|LC_ALL=C awk '{print $3}');
	if [ "$CUSTOM_XRANGE" != "Y" ] ; then
		XMIN=$(head -3 $f|tail -1|LC_ALL=C awk '{print $2}'|\
			perl -e 'while(<>){printf("%g\n",$_);}');
		XMAX=$(head -3 $f|tail -1|LC_ALL=C awk '{print $NF}'|\
			perl -e 'while(<>){printf("%g\n",$_);}');
	fi
	if [ "$CUSTOM_YRANGE" != "Y" ] ; then
		YMIN=$(head -4 $f|tail -1|LC_ALL=C awk '{print $2}'|\
			perl -e 'while(<>){printf("%g\n",$_);}');
		YMAX=$(head -4 $f|tail -1|LC_ALL=C awk '{print $NF}'|\
			perl -e 'while(<>){printf("%g\n",$_);}');
	fi
	VARS=$(head -5 $f|tail -1|sed 's/#//g');
	NC=$(echo $VARS|perl -e 'my @vars,$i=0;while(<>){@vars=split;}
		foreach $v(@vars){if("c" eq $v){print $i,"\n";}++$i;};');
	NPHI=$(echo $VARS|perl -e 'my @vars,$i=0;while(<>){@vars=split;}
		foreach $v(@vars){if("phi" eq $v){print $i,"\n";}++$i;};');
	NPSI=$(echo $VARS|perl -e 'my @vars,$i=0;while(<>){@vars=split;}
		foreach $v(@vars){if("psi" eq $v){print $i,"\n";}++$i;};');
	NGROWTH=$(echo $VARS|perl -e 'my @vars,$i=0;while(<>){@vars=split;}
		foreach $v(@vars){if("phi_growth" eq $v){print $i,"\n";}++$i;};');
	echo "time=$TIME";
	echo "size=${XSIZE}x${YSIZE} [$XMIN:$XMAX]x[$YMIN:$YMAX]";
	echo "vars: C[$NC] PHI[$NPHI] PSI[$NPSI] GROWTH[$NGROWTH]";

	if [ "x$PLOTTER" == "xgnuplot" ] ; then
		gnuplot_plot $f
	elif [ "x$PLOTTER" == "xgri" ] ; then
		gri_plot $f
	else
		echo "unknown plotter: $PLOTTER"
		break
	fi
done

