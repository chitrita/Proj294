# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
C:\Users\David\.spyder2\.temp.py
"""
from __future__ import print_function, division;

import shutil;
import urllib2;
import os;
import sys;
import subprocess as sp;

def read_series_matrix(filename):

    ff = open(filename,'r');
    all_lines = ff.read().splitlines();
    ff.close();
    
        
    exp_dict = dict();
    for line in all_lines:
        split_line = line.split('\t');
        split_line = [sl.replace("\"","") for sl in split_line];
        key = split_line[0];
        key = key.replace('!','');
        if(key == 'Sample_characteristics_ch1'):
           new_key = [x for x in split_line[1:] if len(x) > 0][0].split(':')[0];
           new_values = [''.join(x.split(':')[1:]).strip() for x in split_line[1:]];
           exp_dict[new_key] = new_values;
        for j in range(3,9):
            if(key == ('Sample_supplementary_file_' + str(j))):
                exp_dict[key] = split_line[1:];
            
    return exp_dict;


def select_exp(exp_dict, antibody, time):

    match_antibody = {i for i, x in enumerate(exp_dict['chip antibody']) if antibody.lower() == x.lower()}
    
    corrected_time = str(time) + " min";
           
    match_time = {i for i, x in enumerate(exp_dict['time']) if corrected_time.lower() == x.lower()}
    
    match_sample = match_antibody.intersection(match_time);
    
    if(len(match_sample) == 0):
        print("Error, no matching sample");
        return;
    if(len(match_sample) > 1):
        print("Error, multiple matching samples");
        return;
    
    match_sample = match_sample.pop();
    
    exp_info = dict();
    for key in exp_dict.keys():
        exp_info[key] = exp_dict[key][match_sample];
    
    return exp_info;

def download_sra_files(remote_location, local_location = '', max_recursion = 3, verbose = False):
    """
    Downloads all SRR files in any subdirectory of the remote location
    Max-recursion set to 3 levels (just in case)
    """
  
    def printv(*args):
        if(verbose):
            print(*args);
            sys.stdout.flush();
  
    printv("Reading folder: ", remote_location);  
  
    req = urllib2.Request(remote_location);
    
    response = urllib2.urlopen(req);
    
    the_page = response.read();
    
    entries = the_page.split('\r\n');
    
    #Identify sub folders
    folders = list();
    for entry in entries:
        if(len(entry) == 0):
            continue;
        
        spl_entry = entry.split();
        if(spl_entry[0][0] == 'd'): #if directory flag
            folders.append(spl_entry[-1]);
    
    
    for folder in folders:
        download_sra_files(remote_location + '/' + folder, local_location, max_recursion - 1, verbose);
    
    #Identify SRA files
    files = list();
    for entry in entries:
        if(len(entry) == 0):
            continue;
        
        spl_entry = entry.split();
        if(spl_entry[0][0] == '-' and               #Not a directory
           spl_entry[-1].lower().endswith('.sra')): #Has extension '.sra'
           
           files.append(spl_entry[-1]);
    
    if(len(files) > 0):
        printv("Identified sra files: ");
        for file_name in files:
            printv("   ", file_name);
    
    abs_local_location = os.path.abspath(local_location);
    
    if(not os.path.isdir(abs_local_location)):
        os.makedirs(abs_local_location);
    
    for file_name in files:
    
        printv("Downloading ", file_name);    
    
        file_str = remote_location + '/' + file_name;
    
        req = urllib2.Request(file_str);
        response = urllib2.urlopen(req);
    
        dest_file = open(abs_local_location + os.sep + file_name, 'wb');
        shutil.copyfileobj(response, dest_file)
        dest_file.close();

def sra_to_fastq(sra_filename):
    sra_filename = os.path.abspath(sra_filename);
    sp.check_call(["fastq-dump", sra_filename]);
    return sra_filename.rstrip('.sra')+'.fastq';

#location of the experiment metadata file
exp_file = './Data/GSE36104-GPL15103_series_matrix.txt';

exp_dict = read_series_matrix(exp_file);

antibody = 'IRF4';
time = 0;  #0 30 60 or 120
exp = select_exp(exp_dict, antibody, time);

#Get this from the experiment matrix file
for j in range(3,9):
    key = 'Sample_supplementary_file_' + str(j);
    ftp_str = exp[key];
    if(len(ftp_str) > 0):
        download_sra_files(ftp_str, 'temp', verbose = True);

