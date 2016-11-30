#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
######################
#set dir = $1
export LD_LIBRARY_PATH=/scratch/fkamm/local/lib:$LD_LIBRARY_PATH
export PYTHONPATH="${PYTHONPATH}:/scratch/fkamm/development:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gio:/scratch/fkamm/local/lib/python2.7/site-packages/gtk-2.0/gtk"
source /home/schneider/.bashrc
cd /scratch/schneider/projects/pagerank_refinement/src/experiments/pagerank_additional_test_sets/auto_scripts/

python ../compute_pagerank.py -t /scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/metapsicov_stage2/ -a 0.25 -b 2.5 -l /scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/psicov_lengths -s /scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/psipred/ -o /scratch/schneider/projects/pagerank_refinement/results/pagerank_auto_eval/psicov/metapsicov_stage2/
python ../evaluate_PR.py -t /scratch/schneider/projects/pagerank_refinement/results/pagerank_auto_eval/psicov/metapsicov_stage2/ -a 0.25 -b 2.5 -p /scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/ -o /scratch/schneider/projects/pagerank_refinement/results/pagerank_auto_eval/psicov_metapsicov_stage2_PR.txt -l /scratch/schneider/projects/pagerank_refinement/data/predictor_results/psicov/psicov_lengths
