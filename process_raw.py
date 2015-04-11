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

#Convert input file to FastQ
#returns fastq filename

def ChIPSeq_Pipeline(fastq_file):

    outdirectory = os.path.dirname(fastq_file);
    dumpdirectory = outdirectory + os.sep + '%010x' % random.randrange(16**10)

    sys.path.append('/data/yosef/packages/ChIPSeq/preproc');
    
    #mess with arguments
    
    old_args = sys.argv;
    
    sys.argv = [''];
    sys.argv.append('--fastq_1');
    sys.argv.append(fastq_file);
    
    sys.argv.append('--refgenome_index');
    sys.argv.append('/data/yosef/index_files/mm9/genome/mm9');
    
    sys.argv.append('--out');
    sys.argv.append(dumpdirectory);
    
    chipseqQC = __import__('ChIPSeq-QC');
    
    #Clean up
    sys.argv = old_args;
    target_files = ['aligned_sorted_marked_dumplicates.bam']  #Need to check if this is what we want
    for filename in target_files:    
        shutil.move(dumpdirectory + os.sep + filename, outdirectory);
    
    os.rmdir(dumpdirectory);
    
    return target_files;
    

