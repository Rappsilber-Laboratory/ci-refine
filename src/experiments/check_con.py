import os
import git_tag

file = open("../../data/datasets/ciss_set.txt", "r")
#file = open("../../data/datasets/compiled_sequences_with_seq_dist.txt")

#alphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
#lens = [1.0,1.5,2.0]
#lens = [2.0]
#alphas = [0.4]
lens = [2.0]
alphas = [0.4,0.5]

if os.path.exists("../../results/29-08-14/%s/"%git_tag.get_branch()):
    pass
else:
    os.mkdir("../../results/29-08-14/%s/"%git_tag.get_branch())
for line in file:
    strline = str(line)
    pdb = strline.split()[0][0:5]
    l = strline.split()[1]
    for t in lens:
        for a in alphas:
            print pdb, a, t
            cmd = " ".join(["/scratch/mahmoud/local/bin/python2.7",
                            "/scratch/schneider/projects/rbocon_2.0/src/check_contacts.py",
                            "--restraint_file ../../results/29-08-14/%s/%sRRPAR_%s_%s"%(git_tag.get_branch(),pdb,a,t),
                            #"--restraint_file ../../../rbocon_2.0/results/phy_cm_results/%sRRPCM_1"%(pdb),
                            "--pdb /scratch/schneider/pdb_select_dataset/%s/%s.pdb"%(pdb[0:4],pdb),
                            "--len %s"%l,
                            "--pdb_id %s"%pdb,
                            "--top 2",
                            ">> %s/results_%s_%s.txt"%("../../results/29-08-14/%s/"%git_tag.get_branch(), a, t)])
            os.system( cmd )

