# This module holds some stuff for general coding applications

import sys
import os
#import openbabel

AtomNumList = [(6,4),(7,3),(8,2),(9,1),(15,3),(16,2),(17,1), (35,1),(53,1)]

#usage:get_file_list(string, folderName)
#	Adds all files, containing a given string in name to list
#returns: list with file names
#options: None

#def GetValence(number):
#    valence = 0 
#    for i in AtomNumList:
#        if number = i[0]:
#            valence = i[1]
#    return valence

def get_file_list(string, folderName):
    fileList =[]
    ###Check whether path is valid###
    if os.path.isdir(folderName) == False or os.path.exists(folderName) == False:
        print "No valid path or path doesn`t exists!"
        sys.exit(1)
    ###change dir to folderName###    
    os.chdir(folderName)
    ###create list containing all files in the current folder###
    rawList = os.listdir(os.curdir)
    ###add fileName to list, if it contains the xrayFileKey or the nmrFileKey
    for fileName in rawList:
    	if fileName.find(string) != (-1):
            fileList.append(fileName)
    #for moep in fileList:
    #    print moep
    return fileList
			
			
			
#usage: file2list(file_name)
#	parses a file and adds all lines to a list
#returns: list of strings of the file
#options: None
def file2list(file_name):
	list = []
	file = open(file_name)
	for line in file:
		strline = str(line).strip()
		list.append(strline)
	return list

#def res2mol(file_name, residue):
#   ### store atom numbers
#    file_list = file2list(file_name)
#    file = open("tmp_res", "w")
#    for line in file_list:
#        string = '' 
#        for i in line.split():
#            string += i 
#        #print string, residue.GetName()+str(residue.GetNum())
#        if string.count(residue.GetName()+str(residue.GetNum())) == 1 or string.count(residue.GetName()+'A' + str(residue.GetNum())) == 1 or string.count(residue.GetName()+'B' + str(residue.GetNum())) == 1 or string.count(residue.GetName()+'C' + str(residue.GetNum())) == 1 or string.count(residue.GetName()+'D' + str(residue.GetNum())) == 1:
#            #print "MOEP"
#            file.write(line)
#            file.write('\n')
#    file.close()
#    conv = openbabel.OBConversion()
#    conv.SetInAndOutFormats("pdb", "pdb")
#    mol = openbabel.OBMol()
#    conv.ReadFile(mol, "tmp_res")

#    os.system("rm -rf tmp_res")
#    return mol


