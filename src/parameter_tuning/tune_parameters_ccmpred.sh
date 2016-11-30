#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
######################
#set dir = $1
export LD_LIBRARY_PATH=/scratch/fkamm/local/lib:$LD_LIBRARY_PATH
export PYTHONPATH="${PYTHONPATH}:/scratch/fkamm/development:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gio:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gtk"
source /home/schneider/.bashrc
cd /scratch/schneider/projects/pagerank_refinement/src/
time python tune_pagerank_parameters.py -c ../data/predictor_results/epc-map_test/epc-map_test_length -p /scratch/schneider/projects/pagerank_refinement/data/predictor_results/epc-map_test/ccmpred/ -s /scratch/schneider/projects/pagerank_refinement/data/predictor_results/epc-map_test/psipred/ -m ccmpred
