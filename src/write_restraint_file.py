import InputOutput

restraints = InputOutput.InputOutput.load_restraints("reformatted_HSA_restraints_PR_first_50_perc.txt", seq_sep_min=12)
print restraints
InputOutput.InputOutput.write_contact_file(restraints, "HSA_50_perc_PR_restraints.txt", upper_distance=20)
