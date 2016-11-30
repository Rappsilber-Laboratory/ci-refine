import os
import sys
sys.path.append("../../")
from evaluation_helper import load_pdb_ids, refine_contacts_with_pagerank
from InputOutput import InputOutput

def read_evfold(evfold_file):
    contacts = []
    for line in open(evfold_file):
        strline = str(line).strip().split()
        seq_lower = strline[0]
        seq_upper = strline[2]
        score = strline[5]
        contacts.append((seq_lower, seq_upper, score))
        #sorted(contacts, key=lambda x: x[2])
        #print contacts
        #print len(contacts)
    
    sorted_contacts = sorted(contacts, key=lambda x: x[2], reverse=True)
    #print sorted_contacts
    return sorted_contacts

def main():
    script_folder="/scratch/schneider/projects/pagerank_refinement/src/"
    test_set_folder="/scratch/kstahl/test_sets/rbo_test/"
    output_folder="/scratch/schneider/projects/pagerank_refinement/data/predictor_results/epc-map_test/evfold/"

    protein_data = load_pdb_ids("/scratch/kstahl/test_sets/rbo_test/rbo_test_length")

    os.chdir(script_folder)

    # 1 M 2 E 0.0126265 0

    for pdb_id, length in protein_data:
        contact_file= "".join([test_set_folder,
                             "evfold/",
                             pdb_id,
                             ".evfold"])
       
        contacts = read_evfold(contact_file)
        #break
        #print contacts
        #break
        #os.system("cp %s %s" % (contact_file, output_folder))
        InputOutput.write_contact_file(contacts, "/scratch/schneider/projects/pagerank_refinement/data/predictor_results/epc-map_test/evfold/%s.prediction"%pdb_id)
        
        
main()
