import os
import git_tag

#file = open("../../data/datasets/casp10_len_seq.txt")
file = open("../../data/datasets/ciss_set.txt", "r")
#file = open("../../data/datasets/compiled_sequences_with_seq_dist.txt")

#alphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
#lens = [1.0,1.5,2.0]
#lens = [2.0]
#alphas = [0.4]
lens = [2.0]
alphas = [0.4]
#alphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
#lens = [1.0,1.5,2.0,2.5,3.0]

print "top,acc"
if os.path.exists("../../results/09-09-15/%s/"%git_tag.get_branch()):
    pass
else:
    os.mkdir("../../results/09-09-15/%s/"%git_tag.get_branch())
for line in file:
    strline = str(line)
    pdb = strline.split()[0][0:5]
    l = strline.split()[1]
    for t in lens:
        for a in alphas:
            #print pdb, a, t
            cmd = " ".join(["/scratch/mahmoud/local/bin/python2.7",
                            "/scratch/schneider/projects/rbocon_2.0/src/check_contacts.py",
                            #"--restraint_file /scratch/schneider/projects/rbocon_2.0/results/casp10_RR/RR/%sRR305_1"%(pdb),
                            #"--restraint_file /scratch/schneider/projects/rbocon_2.0/results/casp10_RR/RR/%sRR222_1"%(pdb),
                            #"--restraint_file /scratch/schneider/projects/rbocon_2.0/results/casp10_RR/RR/%sRR113_1"%(pdb),
                            #"--restraint_file /scratch/schneider/projects/rbocon_2.0/results/casp10_RR/RR/%sRR314_1"%(pdb),
                            #"--restraint_file  /scratch/schneider/projects/rbocon_2.0/results/gremlin_pre_casp10/%sRRGRE_1"%(pdb),
                            #"--restraint_file /scratch/schneider/results/graph_kernel_svm/bin_features/psicov_1/gamma_0.001/mi_mid_test/epc_pre_casp/%sRRNone_1"%(pdb),
                            "--restraint_file /scratch/schneider/results/graph_kernel_svm/bin_features/psicov_1/gamma_0.001/mi_mid_test/filtered_results/verify_test/%sRRNone_1"%(pdb),
                            #"--restraint_file ../../results/29-08-14/%s/%sRRPAR_%s_%s"%("merge_sec_struct",pdb,a,t),
                            #"--restraint_file ../../../rbocon_2.0/results/phy_cm_results/%sRRPCM_1"%(pdb),
                            "--pdb /scratch/schneider/pdb_select_dataset/%s/%s.pdb"%(pdb[0:4],pdb),
                            #"--pdb /scratch/schneider/projects/pagerank_refinement/data/pdb/%s.pdb"%(pdb),
                            "--len %s"%l,
                            "--pdb_id %s"%pdb,
                            "--top 2",
                            ">> %s/epc_all_cuts_1.5L_%s_%s.txt"%("../../results/09-09-15/%s/"%git_tag.get_branch(), a, t)])
            os.system( cmd )

