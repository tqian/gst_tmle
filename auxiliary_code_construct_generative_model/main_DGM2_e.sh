#!/bin/bash

cd MISTIE/code # set working directory

FILE_PARA_FOUND="parameter_found.txt"
# indicate finding of the desired parameter value, will be created by parent.R

FILE_LOOP_DONE="DGM2_single_loop_done.txt"
# indicate finish of a batch of parallel jobs, will be created by parent.R

qsub main_DGM2_e_mimicRE_parent.sh # submit the parent job

while [ ! -e "$FILE_PARA_FOUND" ];
do

qsub -t 1:100 main_DGM2_e_mimicRE_parallel.sh # submit parallel jobs

while [ ! -e "$FILE_LOOP_DONE" ];
do
sleep 30 # wait until parallel jobs are finished
done

rm "$FILE_LOOP_DONE" # remove the file for sake of the next loop

done