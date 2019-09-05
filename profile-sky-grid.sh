#!/bin/bash

run_sky_grid() {
	box_size=$1
	image_scale_steps=$2
	/usr/bin/time -l Rscript profile-sky-grid.r $box_size $image_scale_steps > tmp_$1_$2.txt 
        cat tmp_$1_R2.txt \
	  | sed -n 's/\[1\] "\(.*\)"/\1/p; s/.*maximum resident set size  \(.*\)/\1/p; s/.*Elapsed (wall clock) time (h:mm:ss or m:ss): \(.\)/\1/p' \
	  | sed -n '1N; N; s/\n/, /gp'
}

start_experiment() {
	echo
	echo $@
	echo "# im_width, im_height, im_size, im_kbytes, box_size, walltime, maxrss (kb)"
}

# Original image is 356 x 356
# scale=1 is 712x356
# scale=2 is 712x712
# scale=3 is 1424x712
# and so forth...
start_experiment "Profiling image size v/s time for box = 20"
for image_scale_steps in 0 2 4 6 8; do
	run_sky_grid 20 $image_scale_steps
done

start_experiment "Profiling box size v/s memory for image_size = 7.7 Mpix"
for box_size in 10 20 30 40 50 60; do
	run_sky_grid $box_size 6
done

start_experiment "Profiling box size v/s memory for image_size = 30.9 Mpix"
for box_size in 10 20 30 40 50 60; do
	run_sky_grid $box_size 8
done
