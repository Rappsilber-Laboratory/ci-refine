"""Author: Michael Schneider
   Very important notes
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("experiments/")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/structure/")

import helper
import git_tag
from optparse import OptionParser
import numpy
## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-t", type="string", dest="pdb_list", help="Dataset with 5-letter pdb_ids that are used for tuning")
    parser.add_option("-f", type="string", dest="data_path", help="Path with the input contact data")
    parser.add_option("-p", type="string", dest="pdb_path", help="Path with the pdb data")
    parser.add_option("-r", type="string", dest="result_path", help="Path where to write the results to")
    parser.add_option("-s", type="string", dest="psipred_path", help="Path with psipred results")
    options, args = parser.parse_args()
    return options, args 

options, args = add_options(parser)

def pagerank_command(pdb_id,
                     length,
                     top,
                     alpha,
                     hel,
                     pdb_path = options.pdb_path,
                     result_path = options.result_path,
                     sec_struct_path = options.psipred_path,
                     data_path = options.data_path):

    cmd = " ".join(["/scratch/mahmoud/local/bin/python2.7",
                    "/scratch/schneider/projects/pagerank_refinement/src/contact_page_rank.py",
                    "-c %s%sRRGRE_1"%(data_path, pdb_id),
                    "-l %s"%length,
                    "-p %s"%pdb_id,
                    "-f %s%s.pdb"%(pdb_path, pdb_id),
                    "-s %s%s.psipred"%(sec_struct_path, pdb_id),
                    "-t %s"%top,
                    "-a %s"%alpha,
                    "-o %s"%(result_path)])

    return hel.submit_command_job(pdb_id, result_path, cmd, prefix = "j" +  str(alpha) + str(top))

def check_contact_command(pdb_id,
                          length,
                          top,
                          alpha,
                          hel,
                          result_path = options.result_path,
                          pdb_path = options.pdb_path):

    cmd = " ".join(["/scratch/mahmoud/local/bin/python2.7",
                    "/scratch/schneider/projects/rbocon_2.0/src/check_contacts.py",
                    "--restraint_file %s%sRRPAR_%s_%s"%(result_path, pdb_id, alpha, top),
                    "--pdb %s%s.pdb"%(pdb_path, pdb_id),
                    "--len %s"%length,
                    "--pdb_id %s"%pdb_id,
                    "--top 2",
                    "> %s%s_results_%s_%s.txt"%(result_path, pdb_id, alpha, top)])

    return hel.submit_command_job(pdb_id, result_path, cmd, prefix = "j" +  str(alpha) + str(top))

def get_score(result_file):
    score = None
    for line in open(result_file):
        score = float(str(line).split()[8])
    if score == None:
        print score, result_file
        #sys.exit()
        score = 0.5
    return score

def read_pdb_list(pdb_list_file):
    file = open(pdb_list_file)
    pdb_list = []
    for line in file:
        strline = str(line).split()
        pdb_id = strline[0][0:5]
        l = int(strline[1])
        pdb_list.append((pdb_id,l))
    file.close()
    return pdb_list

def main():
    """Generic main function. Executes main functionality of program
    """
    """Some inits"""
    pdb_list = read_pdb_list(options.pdb_list)
    if os.path.exists("%s%s/"%(options.result_path, git_tag.get_branch())):
        pass
    else:
        os.mkdir("%s%s/"%(options.result_path, git_tag.get_branch()))
    res_path = "%s%s/"%(options.result_path, git_tag.get_branch())
    hel = helper.helper()
    jobs = []
    alphas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    tops = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]

    """Run pagerank calculations"""

    for pdb_id, length in pdb_list:
        for alpha in alphas:
            for top in tops:
                job = pagerank_command(pdb_id,
                                       length,
                                       top,
                                       alpha,
                                       hel,
                                       result_path=res_path)
                jobs.append(job)
    hel.wait_for_jobs(jobs)

    jobs = []
    for pdb_id, length in pdb_list:
        for alpha in alphas:
            for top in tops:
                job = check_contact_command(pdb_id,
                                            length,
                                            top,
                                            alpha,
                                            hel,
                                            result_path=res_path)
                jobs.append(job)
    hel.wait_for_jobs(jobs)

    best_score = 0
    best_alpha = None
    best_top = None
    for alpha in alphas:
        for top in tops:
            scores = []
            for pdb_id, _ in pdb_list:
                score = get_score("%s%s_results_%s_%s.txt"%(res_path, pdb_id, alpha, top))
                scores.append(score)
            print alpha, top, numpy.mean(scores)
            if numpy.mean(scores) > best_score:
                best_score = numpy.mean(scores)
                best_alpha = alpha
                best_top = top

    print "Parameter tuning result:", best_score, best_alpha, best_top

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  
