"""
idecoys are intermediate decoys
"""
import random, string
random.seed(26031987)
import pdb
import glob
import os, sys, re, subprocess, shutil
from optparse import OptionParser
import ConfigParser
import src.commands
from time import sleep
import cPickle
import logging
if os.path.exists("template_retrieval/"):
    sys.path.insert(0, "template_retrieval/")
elif os.path.exists("../template_retrieval/"):   
    sys.path.insert(0, "../template_retrieval/")
else:
    sys.path.insert(0,"/scratch/mahmoud/casp_devel/template_retrieval/")
import lib_parser
from itertools import combinations
import time
from email.mime.text import MIMEText
import smtplib
"""Code used in case the decoys in the silent file has number as names
sed -r 's/(^\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)(\S+)$/\1S_\2/' <intermediate_decoys.out > new
sed -r 's/(SCORE:.*\s+)(\S+)$/\1S_\2/' <new >new_intermediate.out
"""

class helper: #Stands for Process INput
    
    abinitio_bin = ""
    cluster_bin = ""
    database_dir = ""
    main_dir = ""
    cparser = None
    parser = None
    converter = None
    def __init__(self):
        self.debug_mode = True

    def openf(self, fname):
        """
        Opens a file for writing, returns the object
        """
        file = open(fname, 'w')
        return file
    
    def wait_for_thread(self, proc):
        """
        Waits for a thread created using subprocess to end
        :arg proc: Objects created by subprocess.Popen()
        :type proc: thread?
        """
        print "Waiting..."
        while proc.poll() is None:
            sleep(5)
            
    def submit_command_job(self, protein, out_folder, command, prefix="j_", queue="casp.q", max_jobs_on_q = 1200, other_options = ""):

        """
        Submit job to the cluster. 
        Input:
        :arg protein: Protein name. The jobname would be prefix+protein
        :arg out_folder:  Folder where to save the shell file and output files
        :arg command: Command to execut on the cluster. Multiple commands can be concatenated by adding \n between each
        :arg prefix: Job name prefix. The final jobname would be prefix+protein
        :arg queue: (all.q or casp.q)
        :arg max_jobs_on_q:  The maximal number of allowed jobs on the queue
        Output:
        the job name. string
        """
        job_name = prefix + protein
        sub_file_name = "%s/%s.sh"%(out_folder, job_name)
        sub_file = open(sub_file_name, 'w')
        sub_file.write("#!/bin/bash\n")
        sub_file.write("#$ -cwd\n")
        sub_file.write("#$ -j y\n")
        sub_file.write("#$ -S /bin/bash\n")
        sub_file.write("######################\n")
        sub_file.write("export LD_LIBRARY_PATH=/scratch/fkamm/local/lib:$LD_LIBRARY_PATH\n")
        sub_file.write("export BAKER_HOME=$BAKER_HOME:/scratch/protein_folding/\n")
        sub_file.write("source /scratch/mahmoud/casp_devel/.bashrc\n")
        sub_file.write("%s\n"%command)
        sub_file.close()
        
        while (self.job_on_q(queue) > max_jobs_on_q):
            sleep(5)
        out_fname = out_folder + "/" + job_name + ".out"
        if os.path.exists(out_fname):
            os.remove(out_fname)
        
        qsub_cmd = " ".join(["qsub",
                             "-q", queue,
                             "-N", job_name,
                             "-o", out_fname,
                             other_options,
                             sub_file_name])
        os.system(qsub_cmd)
        return job_name
        
    def wait_for_jobs(self, jobs):
        """
        Waits for list of jobs to finish. Sleeps until they are finished.
        Input:
        jobs. list(string). List of the jobnames to wait for
        """
        sleep(2)
        print "Waiting ...",
        finished = []
        while len(finished) != len(jobs):
            print ".",
            for job in list(set(jobs)-set(finished)): #iterate over non-finished jobs
                if self.is_job_finished(job):
                    finished.append(job)
            sleep(1)
            finished = list(set(finished)) #Remove duplicate. Just in case
            print len(finished), len(jobs)

    def is_job_finished(self, job_name):
        """Check if the abinitio job is done or not"""
        out = src.commands.getstatusoutput("qstat")[1]
        result = re.findall(job_name[:10] , out) #qstat shows only the first 10 chars of the jobname
        if len(result) > 0:
            return False
        else:
            return True

    def job_on_q(self, queue):
        tmp = src.commands.getstatusoutput('qstat -q ' + queue + ' | wc -l')
        n_jobs = max(int(tmp[1])-2, 0) 
        return n_jobs
    
    def cnfg(self, var, section='S1'):
        """Read config"""
        if self.cparser is None:
            self.cparser= ConfigParser.SafeConfigParser()
            if os.path.exists('config.cfg'):
                self.cparser.read('config.cfg')
            elif os.path.exists('../config.cfg'):
                self.cparser.read('../config.cfg')
            elif os.path.exists('../../config.cfg'):
                self.cparser.read('../../config.cfg')
            elif os.path.exists('../../../config.cfg'):
                self.cparser.read('../../../config.cfg')
            elif os.path.exists('../../../../config.cfg'):
                self.cparser.read('../../../../config.cfg')
            elif os.path.exists('../../../../../config.cfg'):
                self.cparser.read('../../../../../config.cfg')
            else:
                sys.exit("ERROR: config file could not be found. Make sure your are running the script from the main folder or from the compile_test_set folder")
        return self.cparser.get(section, var)
           
    def parse_options(self):
        parser = OptionParser()
        parser.add_option('-p', '--protein', dest='protein', default="1elw", help = 'Query protein.')
        (options, args) = parser.parse_args()
        
        if (options.protein==None):
            parser.print_help()
            sys.exit()
        
        return options

    def extract_silent(self, filename):
        extract_command = self.cluster_bin + \
                          " -in::file::silent " + filename + \
                          " -in::path::database " + self.database_dir
        os.system(extract_command)
        return glob.glob("I_*.pdb")
    
    def download_pdb(self, pdb_id):
        """Download and clean a pdb file)"""
        print "Starting download of %s"%pdb_id
        download_command = self.read_configuration("python_bin") + " " + self.read_configuration("clean_pdb_bin") + " " + pdb_id[:4] + " " + pdb_id[4]
        
        os.system(download_command)
        
        
    def check_protein(self, protein):
        """Check if the protein is:
        - 100<= length(protein) <= 200
        - Has one domain """
        if not (protein[1] == ""): #Multi domain
            return False
        self.download_pdb(protein[0])
        fasta_filename = protein[0][:4]+"_"+protein[0][4]+".fasta"
        pdb_filename = protein[0][:4]+"_"+protein[0][4]+".pdb"
        length = self.get_protein_length(fasta_filename)
        os.system("rm "+ fasta_filename)
        if length >= 100 and length <= 200:
            return True
        else:
            return False
        
    def file2list(self, filename):
        out = []
        try:
            file = open(filename, 'r')
            for line in file:
                if line.strip(): #If length line > 0
                    out.append(line.strip())
        except Exception:
            print "process_input:file2list File %s could not be opened"%(filename)
            sys.exit()
        return out
    
    def mkdir(self, dirname):
        os.system("mkdir -p " + dirname)
        
    def get_protein_length(self, seq_file):
        try:
            with open(seq_file, 'r') as f:
                lines = f.readlines()
                length = 0
                for line in lines[1:]:
                    length = length + len(line.strip())
                return length
        except IOError:
            return 1000
        print "ERROR: Could not open %s"%seq_file
        sys.exit()
        return 1000
        
    def tm_align(self, protein1, protein2):
        """Use Tm align to align two proteins. Returns the alignment length and the TM-score"""
        command = self.read_configuration("tmalign") + " " + \
                protein1 + " " +  protein2
        print command
        
        out = src.commands.getstatusoutput(command)
        #Match TM score
        match = re.search("TM-score=\s+(\S+)\s+.*Chain\_2", out[1])
        if match:
            tm_score = float(match.group(1))
        else:
            print "Error: no alignment possible between %s and %s"%(protein1, protein2)
            print out[1]
            return 0,0.0
        
        #Match alignment length
        match = re.search("Aligned\slength=\s+(\S+),\s+", out[1])
        if match:
            aligned_len = int(match.group(1))
        else:
            print "Error: no alignment possible between %s and %s"%(protein1, protein2)
            print out[1]
            return 0, 0.0
        
        return tm_score, aligned_len
        
    def debug(self, msg):
        if self.debug_mode:
            print "DEBUG: %s"%msg
        return 
    
    def generate_random_str(self, length=5, chars=string.ascii_uppercase):
        """
        http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits
        """
        return ''.join(random.choice(chars) for _ in range(length))

    def setup_logger(self,logger_name, log_file, level=logging.INFO):
        l = logging.getLogger(logger_name)
        formatter = logging.Formatter('%(asctime)s : %(message)s')
        fileHandler = logging.FileHandler(log_file, mode='w')
        fileHandler.setFormatter(formatter)
        streamHandler = logging.StreamHandler()
        streamHandler.setFormatter(formatter)

        l.setLevel(level)
        l.addHandler(fileHandler)
        l.addHandler(streamHandler) 

    def getLogger(self,logger_name):
        return logging.getLogger(logger_name)

    def write_single_fasta_no_line_breaks(self, fasta_seq, fasta_id, output_dir):
        fasta_filename = os.path.join(output_dir, '%s.fasta' %fasta_id)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(">%s\n" %fasta_id)
        # make sure that there is no line break within the fasta_seq
        fasta_seq = fasta_seq.replace("\n","")
        # write fasta_seq with closing line break
        fasta_file.write("%s\n" %fasta_seq)
        fasta_file.close()
        return fasta_filename        
        
    def load_pickle(self, file_name):
        if (os.path.isfile(file_name)):
            try:
                pkl_file = open(file_name, 'rb')
                pkl_content = cPickle.load(pkl_file)
                pkl_file.close()
                return pkl_content
            except Exception, e:
                print 'ERROR: Could not open %s' %(file_name)
                sys.exit()
                return None
        else:
            print 'ERROR: File %s not found\n'%(file_name)
            sys.exit()
            return None   
            
    def dump_pickle(self, file_name, content):
        try:
            with open(file_name, 'w') as f:
                cPickle.dump(content, f, protocol=cPickle.HIGHEST_PROTOCOL)
        except IOError:
            print 'ERROR: Could not open %s' %(file_name)
            sys.exit()

    def align_proteins_according_to_pir(self, pir_fname, target_1, target_2, target_1_seq_fname, target_2_seq_fname, target_1_pdb_fname, target_2_pdb_fname):
        """Aligns two proteins according to their alignment in a pir file. (Note that both proteins
            may be defined as templates in the pir file). Returns the Tm score normalized by the 
            length of the first sequence, second sequence and the length of the alignment..
        """
        pass

    def align_templates_from_pir(self, pir_fname, data_dir, out_fname, do_cache = False):
        """Aligns all pairs of templates in a pir file to each other using
        the alignment in the pir file and the TMalignment.
        Writes a file with the TMscore between each pair of templates
        (normalized by the length of chain 1, normalized by the length of 
        chain 1, normalized by the length of the alignment)"""
        print ">>>> Starting align templates from pir"
        self.parser = lib_parser.Parser()
        self.converter = self.parser.converter
        cache_scores = {}
        cmds = []
        jobs = []
        aln_len_cache = {}
        random.seed(pir_fname + str(time.time())) # initialize RNG
        rand_string = ''.join(random.sample(string.ascii_lowercase, 4))
        tmp_folder = os.path.join(self.cnfg("tmp_dir"), rand_string)
        # print "Get all templates from pir"
        targets = self.parser.get_all_templates_in_pir(pir_fname)
        # print "Parse all pir lines"
        entries = self.parser._parse_pir_lines(pir_fname)
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)
        """ Ready the data """
        combi = list(combinations(targets, 2))
        print len(list(combinations(targets, 2)))
        for i, (elem_1, elem_2) in enumerate(combi):
            #print i, elem_1, elem_2
            # print "Extract the alignment between the two proteins"
            
            rand_string = ''.join(random.sample(string.ascii_lowercase, 4))
            pir_pair_aln = self.parser.extract_custom_pir_pair(pir_fname, 
                                                               elem_1, 
                                                               elem_2, 
                                                               entries)

            # print "Complete the sequences of the two proteins"
            tmp_pir_fname = ".".join([elem_1, elem_2, rand_string, "pir"])
            tmp_pir_fname = os.path.join(tmp_folder, tmp_pir_fname)
            tmp_fasta_fname = ".".join([elem_1, elem_2, rand_string, "fasta"])
            tmp_fasta_fname = os.path.join(tmp_folder, tmp_fasta_fname)
            self.converter.compile_pir_for_one_aln(pir_pair_aln, tmp_pir_fname)
            elem_1_fasta = os.path.join(data_dir, elem_1 + ".fasta")
            elem_2_fasta = os.path.join(data_dir, elem_2 + ".fasta")
            elem_1_pdb = os.path.join(data_dir, elem_1 + ".pdb")
            elem_2_pdb = os.path.join(data_dir, elem_2 + ".pdb")
            self.parser.complete_seq_in_pir(tmp_pir_fname, elem_1, elem_1_fasta, tmp_folder, tmp_pir_fname)
            self.parser.complete_seq_in_pir(tmp_pir_fname, elem_2, elem_2_fasta, tmp_folder, tmp_pir_fname)
            
            # print "Create a fasta file for that alignment"
            shutil.copy(elem_1_pdb, tmp_folder)
            shutil.copy(elem_2_pdb, tmp_folder)
            dir_now = os.getcwd()
            os.chdir(tmp_folder)
            self.converter.pir2fasta(tmp_pir_fname, tmp_fasta_fname)
            os.chdir(dir_now)
            aln_len = self.parser.calc_length_alignment(tmp_fasta_fname)
            aln_len_cache[(elem_1, elem_2)] = aln_len
            cmds.append(self.get_tm_align_cmd(elem_1_pdb, elem_2_pdb, tmp_fasta_fname, aln_len, elem_1, elem_2, tmp_folder)) #-> outfile elem_1_elem_2.tmp
        """ Run the jobs"""
        print len(cmds)
        from math import ceil
        cmd_per_batch = int(ceil(len(cmds)/50.0))
        for i in xrange(0, len(cmds), cmd_per_batch):
            print i, i+cmd_per_batch
            cmd = "\n".join(cmds[i:i+cmd_per_batch])
            rnd_prefix = ''.join(random.sample(string.ascii_lowercase, 6))
            jobs.append(self.submit_command_job( rnd_prefix, # job name
                                                tmp_folder, 
                                                cmd, # the command (concat multiple executions)
                                                prefix="tm", # job name prefix
                                                queue="release.q"))
            # submit_command_job(self, protein, out_folder, command, prefix="j_", queue="casp.q", max_jobs_on_q = 1200, other_options = ""):
            print tmp_folder + "tm_" + rnd_prefix
        self.wait_for_jobs(jobs)
        """ Save the results """
        with open(out_fname, 'w') as out_f:
            for i, (elem_1, elem_2) in enumerate(combi):
                aln_len = 0
                try:
                    aln_len = aln_len_cache[(elem_1, elem_2)]
                    scores = self.parse_tm_scores(elem_1, elem_2, tmp_folder )
                except:
                    scores = [-1, -1, -1]
                out_line = " ".join([elem_1, 
                                     elem_2, 
                                     str(aln_len), 
                                     str(scores[0]), 
                                     str(scores[1]), 
                                     str(scores[2])])
                out_f.write(out_line + "\n")
                print out_line
        shutil.rmtree(tmp_folder)

    def get_tm_align_cmd(self, chain_1_file, chain_2_file, alignment_file, alignment_length, elem_1, elem_2, tmp_folder):
        #print "Calculating TM align scores..."
        # if alignment_length == 0:
        #     return [-1, -1, -1]
        # validate data
        if( not chain_1_file or not os.path.isfile(chain_1_file) ):
            raise IOError("File does not exist: " + str(chain_1_file))
        if( not chain_2_file or not os.path.isfile(chain_2_file) ):
            raise IOError("File does not exist: " + str(chain_2_file))
        if( not alignment_file or not os.path.isfile(alignment_file) ):
            raise IOError("File does not exist: " + str(alignment_file))

        chain_1_file = os.path.abspath(chain_1_file)
        chain_2_file = os.path.abspath(chain_2_file)
        alignment_file = os.path.abspath(alignment_file)

        # validate environment
        if( not self.which('TMalign') ):
            raise RuntimeError('TMalign cannot be found, is it in the PATH?')

        # calculate scores
        cmd = ['TMalign', chain_1_file, chain_2_file, '-I', alignment_file]
        if alignment_length:
            cmd.extend(['-L', str(alignment_length)])
        cmd.extend(['>', os.path.join(tmp_folder, elem_1 + elem_2 + ".tmalign")])
        return " ".join(cmd)

    def parse_tm_scores(self, elem_1, elem_2, tmp_folder ):
        print os.path.join(tmp_folder, elem_1 + elem_2 + ".tmalign")
        with open (os.path.join(tmp_folder, elem_1 + elem_2 + ".tmalign"), "r") as myfile:
            data = myfile.read()
        scores = re.findall("TM-score\s*=\s*([0-9.]+)\s.*", data)
        if not scores:
            print "The output of TMalign could not be parsed:" + data

            raise SyntaxError("The output of TMalign could not be parsed!")
        return scores

    def calc_tm_align_scores(self, chain_1_file, chain_2_file, alignment_file, alignment_length=None ):
        """
        Wraps the TM-align program to perform an alignment and calculate different
        scores.

        Keyword arguments:
        chain_1_file   -- A PDB file
        chain_2_file   -- Another PDB file
        alignment_file -- An FASTA alignment file
        :arg alignment_length: Number of aligned residues
        :type alignment_length: String or int (a number)
        Returns the two TM-scores, normalized by length of chain_1 and chain_2
        Throws exceptions on error.
        Stolen from Philipp's script lib_analysis
        """
        #print "Calculating TM align scores..."
        if alignment_length == 0:
            return [-1, -1, -1]
        # validate data
        if( not chain_1_file or not os.path.isfile(chain_1_file) ):
            raise IOError("File does not exist: " + str(chain_1_file))
        if( not chain_2_file or not os.path.isfile(chain_2_file) ):
            raise IOError("File does not exist: " + str(chain_2_file))
        if( not alignment_file or not os.path.isfile(alignment_file) ):
            raise IOError("File does not exist: " + str(alignment_file))

        chain_1_file = os.path.abspath(chain_1_file)
        chain_2_file = os.path.abspath(chain_2_file)
        alignment_file = os.path.abspath(alignment_file)

        # validate environment
        if( not self.which('TMalign') ):
            raise RuntimeError('TMalign cannot be found, is it in the PATH?')

        # calculate scores
        cmd = ['TMalign', chain_1_file, chain_2_file, '-I', alignment_file]
        if alignment_length:
            cmd.extend(['-L', str(alignment_length)])
        p = subprocess.Popen(cmd , stdout=subprocess.PIPE)
        #print "Executing: %s" % (subprocess.list2cmdline(cmd))
        out, err = p.communicate()

        # handle command execution errors
        if( p.returncode != 0 ):
            raise OSError(p.returncode)

        print out
        scores = re.findall("TM-score\s*=\s*([0-9.]+)\s.*", out)
        if not scores:
            print "The output of TMalign could not be parsed:" + out

            raise SyntaxError("The output of TMalign could not be parsed!")
        return scores

    # http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python/377028#377028
    def which(self, program):

        def is_exe(fpath):
            return os.path.exists(fpath) and os.access(fpath, os.X_OK)

        def ext_candidates(fpath):
            yield fpath
            for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
                yield fpath + ext

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                for candidate in ext_candidates(exe_file):
                    if is_exe(candidate):
                        return candidate

        return None

    def add_to_dict_list(self, my_dict, key, value):
        """Adds a value to a list under key in dict"""
        if key in my_dict:
            my_dict[key].append(value)
        else:
            my_dict[key] = [value]
        return my_dict

    def calculate_alignment_lengths_in_pir(self, pir_fname):
        """Returns a dict containing all templates alignment length to the query"""
        aln_lens = {}
        p = lib_parser.Parser()
        entries = p._parse_pir_lines(pir_fname)
        if entries:
            query_seq = entries[0][2][:-2]
            for i, entry in enumerate(entries[1:]):
                template = entry[0].split(";")[1].strip()
                template_seq = entry[2][:-2]
                aln_len = p.calc_length_alignment_seq(query_seq, template_seq)
                aln_lens[(template, i)] = aln_len
        else:
            print "No entry in pir file %s"%pir_fname
        return aln_lens

    def do_pir_contain_template(self, filename, pdb_name):
        """Check if the pir file have a template in it. Returns true if a pir file which has structural templates exists"""
        if not os.path.exists(filename): 
            return False
        file = open(filename, 'r')
        i = 0
        for line in file:
            if line[0] == ">":
                if (line.strip()).split(";")[1] != pdb_name:
                    #found template other than our target
                    i = i + 1
        file.close()
        return (i >= 1) #True if a template exists in the pir file (more than one sequence in the pir file)
    
    def exec_cmd(self, cmd, log_fname):
        """Runs a command locally and write the log in logfname"""
        print " ".join(cmd)
        if not os.path.exists(os.path.dirname(log_fname)):
            os.system("mkdir -p " + os.path.dirname(log_fname))
        logf = self.openf(log_fname)
        logf.write(" ".join(cmd) + "\n")
        proc = subprocess.Popen(cmd,
                                stdout = logf,
                                stderr = logf)
        self.wait_for_thread(proc)

    def send_email(self, recipients, subject, text, sender = None):
        """
        Sends an email from rbo.deaemon@gmail.com to a list of recipients"""
        if type(recipients) is str:
            recipients = [recipients]

        server = 'smtp.googlemail.com'
        user = 'rbo.daemon'
        password = 'rbocasp2014'
        #In CASP case the acceptance email should be sent to casp-meta AT predictioncenter.org
        if not sender:
            sender =  'rbo.daemon@googlemail.com'
        message = MIMEText(text)
        message['Subject'] = subject
        message['From'] = sender
        message['To'] = recipients[0]
        session = smtplib.SMTP(server)
        session.ehlo()
        session.starttls()
        session.ehlo()
        session.login(user, password)
        session.sendmail(sender, recipients, message.as_string())
        session.quit()


    def add_b_values(self, list_pdbs, target_name, input_folder, output_folder, index = None):
        """add the b_values from the input pdbs to the assembled models because rosetta removes the b_value column when it assembles the domains

            @parameter list_pdbs - list of pdbs containing b_values
            @parameter target_name - target name
            @parameter input_folder - input folder
            @parameter output_folder - output folder
            @parameter index - if index is set index is attached to the file names otherwise the name is exctracted from the file name (only for domains)
        """
        main_dir = self.cnfg("main_dir")
        python_bin = self.cnfg("python_bin", 'S2')
        dic_b_value={}
        for comb in list_pdbs:
            #print comb
     
            with open(os.path.join(input_folder,comb)+'.pdb','r') as rfile:
                old_id = None
                for line in rfile:
                    #line_list=filter(None,line.split(' '))
                    
                    if line[:4]=='ATOM':
                        
                        residue_id=line[22:26]
                        if residue_id != old_id:
                            if index != None:
                                dic_b_value.setdefault(int(index), []).append(line[60:66].strip())
                            else:
                                dic_b_value.setdefault(int(comb.split('_')[-2][1:]),[]).append(line[60:66].strip())
                        old_id = residue_id
                  
        list_b_values=[]
        for dom in sorted(dic_b_value.keys()):#make sure that the order of the domains is from the first to the last
            for b_value in dic_b_value[dom]: 
                list_b_values.append(b_value)
        #print dic_b_value
        txt=[]
        for index,b_value in enumerate(list_b_values):
            txt.append(str(index+1))
            txt.append('\t')
            txt.append(b_value)
            txt.append('\n')
                   
            with open(os.path.join(output_folder,'%s_b_values.txt'%target_name),'w') as fopen:
               fopen.write(''.join(txt))
            
            add_bvalue_cmd = " ".join([python_bin,
                                     os.path.join(main_dir,'tools/pdb_bfactor.py'),
                                     '-i', os.path.join(output_folder,'%s.pdb'%target_name),
                                     '-d', os.path.join(output_folder,'%s_b_values.txt'%target_name),
                                     '-o', os.path.join(output_folder,target_name+'_tmp_bvalue.pdb')])
            #print add_bvalue_cmd
            os.system(add_bvalue_cmd)
            shutil.copy(os.path.join(output_folder,target_name+'_tmp_bvalue.pdb'),os.path.join(output_folder,'%s.pdb'%target_name))

