#!/bin/bash 
#PJM --rsc-list "node=1"  
#PJM --rsc-list "rscunit=rscunit_ft01"
#PJM --rsc-list "rscgrp=eap-small"  
#PJM --rsc-list "elapse=60:00" 
#PJM --mpi "proc=4"
#PJM -S   

module load lang

export OMP_NUM_THREADS=12

fapp -C -Hevent=pa1 -d ./rep1 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa2 -d ./rep2 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa3 -d ./rep3 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa4 -d ./rep4 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa5 -d ./rep5 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa6 -d ./rep6 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa7 -d ./rep7 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa8 -d ./rep8 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa9 -d ./rep9 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa10 -d ./rep10 -L 4 mpiexec ./monet_k
fapp -C -Hevent=pa11 -d ./rep11 -L 4 mpiexec ./monet_k
