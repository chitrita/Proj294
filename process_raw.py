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

def ChIPSeq_Pipeline_Single(fastq_file, dumpdirectory=''):

    if(not os.path.isfile(fastq_file)):
        print("Error: File does not exist:", fastq_file);
        return;

    fastq_file = os.path.abspath(fastq_file);
    outdirectory = os.path.dirname(fastq_file);


    if(dumpdirectory == ''):
        bname = os.path.basename(fastq_file).rstrip(".fastq");
        dumpdirectory = outdirectory + os.sep + bname + '_ChIPSeq_' + '%010x' % random.randrange(16**10);

    if(not os.path.isdir(dumpdirectory)):
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

    old_dir = os.getcwd();
    os.chdir(os.path.dirname(CS_PIPELINE_LOCATION));
    sys.path.insert(0, os.path.dirname(CS_PIPELINE_LOCATION));
    chipseqQC = execfile(CS_PIPELINE_LOCATION, CS_globals, CS_locals);

    #Clean up
    os.chdir(old_dir);
    sys.argv = old_args;
    sys.path.remove(os.path.dirname(CS_PIPELINE_LOCATION));

def ChIPSeq_Pipeline_Paired(fastq_file1, fastq_file2, dumpdirectory=''):


    fastq_file1 = os.path.abspath(fastq_file1);
    fastq_file2 = os.path.abspath(fastq_file2);

    if(not os.path.isfile(fastq_file1)):
        print("Error: File does not exist:", fastq_file1);
        return;

    if(not os.path.isfile(fastq_file2)):
        print("Error: File does not exist:", fastq_file2);
        return;

    outdirectory = os.path.dirname(fastq_file1);

    if(dumpdirectory == ''):
        bname = os.path.basename(fastq_file1).rstrip("_1.fastq");
        dumpdirectory = outdirectory + os.sep + bname + '_ChIPSeq_' + '%010x' % random.randrange(16**10);

    if(not os.path.isdir(dumpdirectory)):
        os.makedirs(dumpdirectory);


    #mess with arguments

    old_args = sys.argv;

    sys.argv = [''];
    sys.argv.append('--fastq_1');
    sys.argv.append(fastq_file1);


    sys.argv.append('--fastq_2');
    sys.argv.append(fastq_file2);

    sys.argv.append('--refgenome_index');
    sys.argv.append('/data/yosef/index_files/mm9/genome/mm9');

    sys.argv.append('--paired');

    sys.argv.append('--out');
    sys.argv.append(dumpdirectory);

    old_dir = os.getcwd();
    os.chdir(os.path.dirname(CS_PIPELINE_LOCATION));
    sys.path.insert(0, os.path.dirname(CS_PIPELINE_LOCATION));

    chipseqQC = execfile(CS_PIPELINE_LOCATION, CS_globals, CS_locals);


    #Clean up
    os.chdir(old_dir);
    sys.argv = old_args;
    sys.path.remove(os.path.dirname(CS_PIPELINE_LOCATION));

def remove_intermediates(dumpdirectory = ''):
    if dumpdirectory == '':
        dumpdirectory = os.getcwd();

    files_to_remove = ["initial_bowtie_alignment.bam",
                       "initial_bowtie_alignment.sam"];

    for filename in files_to_remove:
        os.remove(dumpdirectory + os.sep + filename);

if (__name__ == "__main__"):
    if(len(sys.argv) == 3):
        fastq_file = sys.argv[1];
        dumpdirectory = sys.argv[2];
        ChIPSeq_Pipeline_Single(fastq_file, dumpdirectory);
        os.remove(fastq_file);
        remove_intermediates(dumpdirectory);
    if(len(sys.argv) == 4):
        fastq_file1 = sys.argv[1];
        fastq_file2 = sys.argv[2];
        dumpdirectory = sys.argv[3];
        ChIPSeq_Pipeline_Paired(fastq_file1, fastq_file2, dumpdirectory);
        os.remove(fastq_file1);
        os.remove(fastq_file2);
        remove_intermediates(dumpdirectory);
