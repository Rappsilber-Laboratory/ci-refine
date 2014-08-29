import os 

file = open("ciss_set.txt")

#alphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
#lens = [1,1.5,2]
lens = [5.0]
alphas = [0.4]
for line in file:
    strline = str(line).strip().split()
    pdb = strline[0]
    l = strline[1]
    for t in lens:
        for a in alphas:
            print strline
            os.system("/scratch/fkamm/local/bin/python2.7  contact_page_rank.py -c /scratch/schneider/results/graph_kernel_svm/bin_features/psicov_1/gamma_0.001/mi_mid_test/filtered_results/verify_test/%sRRNone_1 -l %s -p %s -f /scratch/schneider/pdb_select_dataset/%s/%s.pdb -s /scratch/schneider/projects/rbocon_2.0/data/rbo_test/sec_struct/%s.psipred -t %s -a %s"%(pdb,l,pdb,pdb[0:4],pdb,pdb,t,a))
    #break
file.close()
