#!/bin/bash

run_profound() {
	box_size=$1
	image_scale_steps=$2
	/usr/bin/time -l Rscript profile.r $box_size $image_scale_steps $what > tmp.txt 2>tmp1.txt
        echo -n `cat tmp.txt`| sed -e's|"||g' > tmpa.txt
        echo -n `grep "real" tmp1.txt` | sed -e"s|  | |g" | cut  -f 1 -d" " > tmpb.txt
        echo -n `grep "maximum resident set size" tmp1.txt`| sed -e"s|maximum resident set size||" > tmpc.txt
        echo `cut -c4-100 tmpa.txt`,`cat tmpb.txt`,`cat tmpc.txt`

}

start_experiment() {
	echo 1>&2
	echo $@ 1>&2
	echo "# im_width, im_height, im_size, im_kbytes, box_size, walltime, maxrss"
}

if [ $# -lt 1 ]; then
	cat <<EOF
Usage: $0 <what>

<what> is one of "profound-original", "profound-adacs", "skygrid-original" and "skygrid-adacs"
EOF
	exit 1
fi
what=$1

# Original image is 356 x 356
# scale=1 is 712x356
# scale=2 is 712x712
# scale=3 is 1424x712
# and so forth...
rm -f bmask-$what.csv
start_experiment "Profiling image size v/s time for box = 20" >> bmask-$what.csv
for image_scale_steps in 0 2 4 6 8; do
	run_profound 20 $image_scale_steps >> bmask-$what.csv
done

start_experiment "Profiling box size v/s memory for image_size = 7.7 Mpix"
for box_size in 10 20 30 40 50 60; do
	run_profound $box_size 6 >> bmask-$what.csv
done

start_experiment "Profiling box size v/s memory for image_size = 30.9 Mpix"
for box_size in 10 20 30 40 50 60; do
	run_profound $box_size 8 >> bmask-$what.csv
done
