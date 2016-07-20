import os
import sys

methods = ["epc-map", "metapsicov_stage2", "evfold"]
datasets = ["psicov", "metapsicov_test", "d329", "svmcon"]

# set up data folders
root_folder = "/scratch/schneider/projects/pagerank_refinement/results/pagerank_auto_eval"
for d in datasets:
    for m in methods:
        experiment_folder = "/".join([root_folder,           
                                      d,     
                                      m])
        os.system("mkdir -p %s" % experiment_folder)


os.system("qsub -o /scratch/schneider/run_tmp/ epc-map_psicov.sh")
os.system("qsub -o /scratch/schneider/run_tmp/ evfold_psicov.sh")
os.system("qsub -o /scratch/schneider/run_tmp/ metapsicov_stage2_psicov.sh")

os.system("qsub -o /scratch/schneider/run_tmp/ epc-map_svmcon.sh")       
os.system("qsub -o /scratch/schneider/run_tmp/ evfold_svmcon.sh")
os.system("qsub -o /scratch/schneider/run_tmp/ metapsicov_stage2_svmcon.sh")        

os.system("qsub -o /scratch/schneider/run_tmp/ epc-map_d329.sh")       
os.system("qsub -o /scratch/schneider/run_tmp/ evfold_d329.sh")
os.system("qsub -o /scratch/schneider/run_tmp/ metapsicov_stage2_d329.sh")  

os.system("qsub -o /scratch/schneider/run_tmp/ epc-map_metapsicov_test.sh")       
os.system("qsub -o /scratch/schneider/run_tmp/ evfold_metapsicov_test.sh")
os.system("qsub -o /scratch/schneider/run_tmp/ metapsicov_stage2_metapsicov_test.sh")  #
