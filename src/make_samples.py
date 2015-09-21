import os 
import sys 

#samples = ["0.25", "0.5", "0.75", "1.0", "1.25", "1.5"]
samples = ["1.75"]
for s in samples:
    for i in xrange(0, 10): 
        max_links = int(578*float(s))
        os.system("/scratch/mahmoud/local/bin/python2.7 xl_page_rank_backup_sort_of_working.py -e /scratch/schneider/projects/pagerank_refinement/data/xl_data/20PercentFDR_xiFDR0.csv -o -28 -n HSA_%sL_%s -m %s"%(s,i, max_links))
        os.system("/scratch/mahmoud/local/bin/python2.7 /scratch/schneider/Software/epc-map/rbocon/release/analysis/write_checked_distances.py --pdb /scratch/schneider/projects/xl_project/data/hsa/1ao6/1ao6A.pdb --pdb_id 1ao6A --len %s --dom_ass 1-578 --atom CA --restraint_file HSA_%sL_%s_PSM.txt"%(int(max_links/2.0), s, i)) 
        os.system("/scratch/mahmoud/local/bin/python2.7 /scratch/schneider/Software/epc-map/rbocon/release/analysis/write_checked_distances.py --pdb /scratch/schneider/projects/xl_project/data/hsa/1ao6/1ao6A.pdb --pdb_id 1ao6A --len %s --dom_ass 1-578 --atom CA --restraint_file HSA_%sL_%s_PR.txt"%(int(max_links/2.0), s, i))

