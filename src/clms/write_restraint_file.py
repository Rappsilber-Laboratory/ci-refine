import src.InputOutput

# 1aoaA
restraints = src.InputOutput.InputOutput.load_restraints("../data/input_xl_refined_data_prediction/reformatted_20_FDR_HSA_restraints_PR_first_50_perc.txt", seq_sep_min=12, pos_min=1, pos_max=197, offset=0)
src.InputOutput.InputOutput.write_restraint_file_hack(restraints, "../data/input_xl_refined_data_prediction/1aoaA_20_FDR_50_perc_PR_restraints.txt", upper_distance=20)

# 1aob
restraints = src.InputOutput.InputOutput.load_restraints("../data/input_xl_refined_data_prediction/reformatted_20_FDR_HSA_restraints_PR_first_50_perc.txt", seq_sep_min=12, pos_min=198, pos_max=386, offset=-197)
src.InputOutput.InputOutput.write_restraint_file_hack(restraints, "../data/input_xl_refined_data_prediction/1aobA_20_FDR_50_perc_PR_restraints.txt", upper_distance=20)

# 1aocA
restraints = src.InputOutput.InputOutput.load_restraints("../data/input_xl_refined_data_prediction/reformatted_20_FDR_HSA_restraints_PR_first_50_perc.txt", seq_sep_min=12, pos_min=387, pos_max=578, offset=-386)
src.InputOutput.InputOutput.write_restraint_file_hack(restraints, "../data/input_xl_refined_data_prediction/1aocA_20_FDR_50_perc_PR_restraints.txt", upper_distance=20)
print restraints


