import os 

file = open("ciss_set.txt", "r")

#alphas = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
#lens = [1.0,1.5,2.0]
lens = [5.0]
alphas = [0.4]
for line in file:
    strline = str(line)
    pdb = strline.split()[0][0:5]
    l = strline.split()[1]
    for t in lens:
        for a in alphas:

            os.system("/scratch/fkamm/local/bin/python2.7 /scratch/schneider/projects/rbocon_2.0/src/check_contacts.py --restraint_file results/%sRRPAR_%s_%s --pdb /scratch/schneider/pdb_select_dataset/%s/%s.pdb --len %s --pdb_id %s --top %s >> results_%s_%s.txt"%(pdb,a,t, pdb[0:4],pdb, l, pdb, 2,a,t))
            #os.system("/scratch/fkamm/local/bin/python2.7 /scratch/schneider/projects/rbocon_2.0/src/check_contacts.py --restraint_file %sRRNone_1 --pdb /scratch/schneider/pdb_select_dataset/%s/%s.pdb --len %s --pdb_id %s --top %s >> epc_results.txt"%(pdb, pdb[0:4],pdb, l, pdb, 2))
