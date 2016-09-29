#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
######################
#set dir = $1
export LD_LIBRARY_PATH=/scratch/fkamm/local/lib:$LD_LIBRARY_PATH
export PYTHONPATH="${PYTHONPATH}:/scratch/fkamm/development:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gio:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gtk"
source /home/schneider/.bashrc
cd /scratch/schneider/projects/pagerank_refinement/src/experiments/pagerank_additional_test_sets/auto_scripts_original/

#python ../compute_pagerank.py -t /scratch/schneider/projects/pagerank_refinement/data/predictor_results/metapsicov_test/metapsicov_stage2/ -a 0.25 -b 2.5 -l /scratch/schneider/projects/pagerank_refinement/data/predictor_results/metapsicov_test/metapsicov_test_lengths -s /scratch/schneider/projects/pagerank_refinement/data/predictor_results/metapsicov_test/psipred/ -o /scratch/schneider/projects/pagerank_refinement/metapsicov_test/metapsicov_stage2/
python ../evaluate_predictions.py -t /scratch/schneider/projects/pagerank_refinement/data/predictor_results/metapsicov_test/metapsicov_stage2/ -a 0.25 -b 2.5 -p /scratch/schneider/projects/pagerank_refinement/data/predictor_results/data/predictor_results/metapsicov_test/ -o /scratch/schneider/projects/pagerank_refinement/metapsicov_test_metapsicov_stage2_PR.txt -l /scratch/schneider/projects/pagerank_refinement/data/predictor_results/metapsicov_test/metapsicov_test_lengths
