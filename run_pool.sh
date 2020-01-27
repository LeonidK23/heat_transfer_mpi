#!/bin/bash
for ghost_lines in `seq 2 2 40`;
do
  echo "Ghost lines size: $ghost_lines"
  mpirun -np 17 -hostfile Hosts.txt ./stencil 1024 $ghost_lines 400 0.4 2 1
done
