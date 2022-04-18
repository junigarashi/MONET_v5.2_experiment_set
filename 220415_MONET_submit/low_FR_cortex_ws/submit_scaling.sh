#!/bin/bash

for dirlist2 in `ls -l | awk '$1 ~ /d/ {print $9 }' `
do
    echo $dirlist2
    cd $dirlist2
    python make_process_configuration.py
    pjsub run_monet_fugaku.sh
    cd ..
done



