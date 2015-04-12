# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:06:31 2015

@author: David
"""

import sys;
import os;
import subprocess as sp;
import random;
import shutil;

CS_PIPELINE_LOCATION = '/data/yosef/packages/ChIPSeq/preproc/ChIPSeq-QC.py';
CS_globals = {'__file__':CS_PIPELINE_LOCATION};
CS_locals = {};


#Convert input file to FastQ
#returns fastq filename

def ChIPSeq_Pipeline(fastq_file):

    if(not os.path.isfile(fastq_file)):
        print("Error: File does not exist:", fastq_file);
        return;

    outdirectory = os.path.dirname(fastq_file);
    dumpdirectory = outdirectory + os.sep + '%010x' % random.randrange(16**10)
    os.makedirs(dumpdirectory);
    
    
    #mess with arguments
    
    old_args = sys.argv;
    
    sys.argv = [''];
    sys.argv.append('--fastq_1');
    sys.argv.append(fastq_file);
    
    sys.argv.append('--refgenome_index');
    sys.argv.append('/data/yosef/index_files/mm9/genome/mm9');
    
    sys.argv.append('--out');
    sys.argv.append(dumpdirectory);
    
    chipseqQC = execfile(CS_PIPELINE_LOCATION, CS_globals, CS_locals);
    
    #Clean up
    sys.argv = old_args;
    target_file = 'aligned_sorted_marked_duplicates.bam';  #Need to check if this is what we want   
    new_file_name = outdirectory + os.sep + os.path.basename(fastq_file).rstrip('.fastq') + '.bam';
    shutil.move(dumpdirectory + os.sep + target_file, new_file_name);
    
    shutil.rmtree(dumpdirectory);
 
    return new_file_name;
    

