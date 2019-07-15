#!/bin/bash

cd MISTIE/code_determine_nmax/exW # set working directory

FILE_PARA_FOUND="nmax_found.txt"
# indicate finding of the desired parameter value, will be created by parent.R

FILE_LOOP_DONE="nmax_single_loop_done.txt"
# indicate finish of a batch of parallel jobs, will be created by parent.R

qsub determine_nmax_parent.sh # submit the parent job

while [ ! -e "$FILE_PARA_FOUND" ];
do

qsub -t 1:500 determine_nmax_parallel.sh # submit parallel jobs

while [ ! -e "$FILE_LOOP_DONE" ];
do
sleep 30 # wait until parallel jobs are finished
done

rm "$FILE_LOOP_DONE" # remove the file for sake of the next loop

done