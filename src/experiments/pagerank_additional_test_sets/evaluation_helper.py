import os 

def load_pdb_ids(pdb_id_file):
    protein_data = []
    for line in open(pdb_id_file):
        pdb_id = str(line).strip().split()[0]
        length = str(line).strip().split()[1]
        protein_data.append((pdb_id, length))
    return protein_data


def evaluate_contact_file(pdb_file, pdb_id, contact_file, length, dom_ass):
    cmd = " ".join(["python",
                    "check_contacts.py",
                    "--pdb", pdb_file,
                    "--pdb_id", pdb_id,
                    "--restraint_file", contact_file,
                    "--len", length,
                    "--dom_ass", dom_ass])
    os.system(cmd)


def refine_contacts_with_pagerank(pdb_id, contact_file, length, psipred_file, top_ratio, output_folder):
    cmd = " ".join(["python",
                    "contact_page_rank.py",
                    "-p", pdb_id,
                    "-c", contact_file,
                    "-l", str(length),
                    "-s", psipred_file,
                    "-t", top_ratio,
                    "-o", output_folder])
    print cmd
    os.system(cmd)

