#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#export LD_LIBRARY_PATH=/scratch/fkamm/local/lib:$LD_LIBRARY_PATH
#export PYTHONPATH="${PYTHONPATH}:/scratch/fkamm/development:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gio:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gtk"
cd /scratch/schneider/projects/pagerank_refinement/src/
source /scratch/schneider/projects/casp_infra/src/casp/.bashrc
#which python
/scratch/mahmoud/local/bin/python2.7 contact_patterns.py  -p ../data/datasets/h_1.txt
#/scratch/mahmoud/local/bin/python2.7 contact_patterns.py  -p ../data/datasets/epc_map_trainset_before_casp10.txt
