import os 
import subprocess


def get_branch():
    p = subprocess.Popen(["/opt/rocks/bin/git rev-parse --abbrev-ref HEAD"],stdout=subprocess.PIPE,shell=True)
    return p.communicate()[0].strip()
    
#print get_branch()
