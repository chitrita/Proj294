from __future__ import division, print_function;

import meta_data;
import file_operations;
import os;
import subprocess as sp;
import re;

LOCATION = os.path.dirname(os.path.abspath(__file__));

def create_directories():
    for exp in meta_data.exp_list:
        directory = meta_data.directory_from_exp(exp);
        os.makedirs(directory);

def download_raw_reads():
    for i,exp in enumerate(meta_data.exp_list):
        print(i, exp["chip antibody"], exp["time"]);
        file_operations.download_sra_for_sample(exp);

def convert_sra_to_fastq():
    for i, exp in enumerate(meta_data.exp_list):
        print(i, exp["chip antibody"], exp["time"]);
        file_operations.reads_to_fastq_for_sample(exp);

def map_reads_for_exp(exp, location = ''):
    """
    Finds the fastq files for an experiment.
    Generates a run.sh file for the queue that calls process_raw
    Submits the run.sh to the queue
    """
    if(location == ''):
        exp_dir = meta_data.directory_raw_reads(exp);
    else:
        exp_dir = location;

    all_files = os.listdir(exp_dir);

    isPaired = False;

    for filename in all_files:
        if(filename.endswith('_1.fastq')):
            isPaired = True;
            fastq1 = exp_dir + os.sep + filename;
        elif(filename.endswith('_2.fastq')):
            isPaired = True;
            fastq2 = exp_dir + os.sep + filename;
        elif(filename.endswith('.fastq')):
             isPaired = False;
             fastq = exp_dir + os.sep + filename;

    with open(LOCATION + os.sep + 'run_CS_pipeline.sh', 'r') as fin:
        run_contents = fin.read();

    if(isPaired):
        run_contents = run_contents.replace("<#FASTQ_FILE1>", fastq1);
        run_contents = run_contents.replace("<#FASTQ_FILE2>", fastq2);
    else:
        run_contents = run_contents.replace("<#FASTQ_FILE1>", fastq);
        run_contents = run_contents.replace("<#FASTQ_FILE2> ", "");

    run_contents = run_contents.replace("<#DUMPDIRECTORY>", exp_dir);
    run_contents = run_contents.replace("<#LOCATION>", LOCATION);

    run_filename = exp_dir + os.sep + "run_CS_pipeline.sh";
    with open(run_filename, 'w') as fout:
        fout.write(run_contents);

    sp.check_call(["qsub", run_filename]);

def manorm_for_exp_pair(exp1, exp2):
   command = "python";
   command = command + " " + LOCATION + os.sep + "manorm.py";
   command = command + " " + meta_data.directory_from_exp(exp1);
   command = command + " " + meta_data.directory_from_exp(exp2);

   if(not os.path.isdir(meta_data.directory_manorm(exp2))): 
      os.mkdir(meta_data.directory_manorm(exp2));

   file_operations.submit_queue(meta_data.directory_manorm(exp2), command); 

def deseq_for_exp(exp):

    time = re.search("\d+", exp["time"]).group();
    antibody = exp["chip antibody"];

    command = "python";
    command = command + " " + LOCATION + os.sep + "deseq2.py";
    command = command + " " + antibody; 
    command = command + " " + time; 

    if(not os.path.isdir(meta_data.directory_deseq(exp))):
        os.mkdir(meta_data.directory_deseq(exp));

    file_operations.submit_queue(meta_data.directory_deseq(exp), command);


def foldchange_for_exp_pair(exp1, exp2):
    command = "python";
    command = command + " " + LOCATION + os.sep + "foldchange.py";
    command = command + " " + meta_data.directory_from_exp(exp1);
    command = command + " " + meta_data.directory_from_exp(exp2);

    if(not os.path.isdir(meta_data.directory_foldchange(exp2))):
        os.mkdir(meta_data.directory_foldchange(exp2));

    file_operations.submit_queue(meta_data.directory_foldchange(exp2), command);

def macs2_for_exp(exp):
    command = "python";
    command = command + " " + LOCATION + os.sep + "macs2.py";
    command = command + " " + meta_data.directory_from_exp(exp);

    if(not os.path.isdir(meta_data.directory_peaks(exp))):
        os.mkdir(meta_data.directory_peaks(exp));

    file_operations.submit_queue(meta_data.directory_peaks(exp), command);

