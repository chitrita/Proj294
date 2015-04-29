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
import fileinput;
import meta_data;


def isPairedSRA(filename):
    filename = os.path.abspath(filename);
    try:
        with open(os.devnull, "w") as fnull:
            contents = sp.check_output(["fastq-dump","-X","1","-Z","--split-spot", filename], stderr=fnull);
    except sp.CalledProcessError, e:
        raise Exception("Error running fastq-dump on",filename);

    if(contents.count("\n") == 4):
        return False;
    elif(contents.count("\n") == 8):
        return True;
    else:
        raise Exception("Unexpected output from fast-dump on ", filename);



def download_sra_files(remote_location, local_location = '', max_recursion = 3, verbose = False):
    """
    Downloads all SRR files in any subdirectory of the remote location
    Max-recursion set to 3 levels (just in case)
    """

    downloaded_files = list();

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
        dl_files = download_sra_files(remote_location + '/' + folder, local_location, max_recursion - 1, verbose);
        downloaded_files.extend(dl_files);

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

        dest_file_name = abs_local_location + os.sep + file_name;
        dest_file = open(dest_file_name, 'wb');
        shutil.copyfileobj(response, dest_file)
        dest_file.close();
        downloaded_files.append(dest_file_name);

    return downloaded_files;

def sra_single_to_fastq(sra_filename):
    sra_filename = os.path.abspath(sra_filename);
    try:
        sp.check_call(["fastq-dump",'--outdir',os.path.dirname(sra_filename), sra_filename]);
        return sra_filename.rstrip('.sra')+'.fastq';
    except sp.CalledProcessError, e:
        print("Error reading file:", sra_filename);
        return None;


def sra_paired_to_fastq(sra_filename):
    sra_filename = os.path.abspath(sra_filename);
    try:
        sp.check_call(["fastq-dump", "--split-files", "--outdir", os.path.dirname(sra_filename), sra_filename]);
        base = sra_filename.rstrip('.sra');
        return (base+'_1.fastq', base+'_2.fastq');
    except sp.CalledProcessError, r:
        print("Error reading file:", sra_filename);
        return None;

def bam_to_bed(bam_filename):
    bam_filename = os.path.abspath(bam_filename);
    bed_filename = bam_filename.rstrip('.bam') + '.bed';
    with open(bed_filename, 'w') as fout:
        sp.check_call(['bedtools', 'bamtobed', '-i', bam_filename], stdout=fout)

    return bed_filename;

def download_sra_for_sample(exp, local_location=''):
    if(local_location == ''):
        local_location = meta_data.directory_raw_reads(exp);

    sra_files = list();

    for j in range(1,10):
        key = 'Sample_supplementary_file_' + str(j);
        ftp_str = exp[key];
        if(len(ftp_str) > 0 and "sra" in ftp_str):
            sra_files.extend(download_sra_files(ftp_str, local_location, verbose = True));

    return sra_files;

def reads_to_fastq_for_sample(exp, loc = ''):
    if(loc == ''):
        loc = meta_data.directory_raw_reads(exp);
    antibody = exp["chip antibody"];
    time = exp["time"].replace(" ", "");
    sra_files = [loc + os.sep + filename for filename in os.listdir(loc) if filename.lower().endswith('.sra')];

    if(isPairedSRA(sra_files[0])):  #Assumes if the first file is paired, then they all are
        fastq_1_files = list();
        fastq_2_files = list();
        for filename in sra_files:
            f1,f2 = sra_paired_to_fastq(filename);
            if f1 is not None: fastq_1_files.append(f1);
            if f2 is not None: fastq_2_files.append(f2);

        out_fastq1_name = str(antibody) + "_" + str(time) + "_1.fastq";
        out_fastq1_name = loc + os.sep + out_fastq1_name;


        out_fastq2_name = str(antibody) + "_" + str(time) + "_2.fastq";
        out_fastq2_name = loc + os.sep + out_fastq2_name;

        with open(out_fastq1_name, 'w') as fout:
            for line in fileinput.input(fastq_1_files):
                fout.write(line);

        with open(out_fastq2_name, 'w') as fout:
            for line in fileinput.input(fastq_2_files):
                fout.write(line);

        for filename in fastq_1_files:
            os.remove(filename);
        for filename in fastq_2_files:
            os.remove(filename);
        for filename in sra_files:
            os.remove(filename);

        return [out_fastq1_name, out_fastq2_name];

    else:  #Unpaired SRA files

        fastq_files = list();
        for filename in sra_files:
            f = sra_single_to_fastq(filename);
            if f is not None: fastq_files.append(f);

        out_fastq_name = str(antibody) + "_" + str(time) + ".fastq";
        out_fastq_name = loc + os.sep + out_fastq_name;
        with open(out_fastq_name, 'w') as fout:
            for line in fileinput.input(fastq_files):
                fout.write(line)

        #Clean up partial fastq files
        for filename in fastq_files:
            os.remove(filename);
        for filename in sra_files:
            os.remove(filename);

        return out_fastq_name;

def download_peaks_for_sample(exp, outputdir = ''):
    if(outputdir == ''):
        outputdir = meta_data.directory_peaks(exp);

    abs_outputdir = os.path.abspath(outputdir);

    if(not os.path.isdir(abs_outputdir)):
        os.mkdir(abs_outputdir);

    file_str = exp['Sample_supplementary_file_2'];
    file_name = os.path.basename(file_str);

    print("Downloading ", file_str);

    req = urllib2.Request(file_str);
    response = urllib2.urlopen(req);

    dest_file_name = abs_outputdir + os.sep + file_name;
    dest_file = open(dest_file_name, 'wb');
    shutil.copyfileobj(response, dest_file)
    dest_file.close();

    #Unzip if necessary
    if(dest_file_name.endswith('.gz')):
        sp.check_call(['gzip','-d',dest_file_name]);
        dest_file_name = dest_file_name.rstrip('.gz');

    return dest_file_name;

def bed_merge(output_file, *inputfiles):
    """
    Runs bedtools merge command on all input files.
    Sends output to output_file
    Returns output_File

    If output_file contains a path (at least one path seperator charactrer)
       then output to that location.  Otherwise, put in same location as the
       first input file.
    """
    working_dir = os.path.dirname(inputfiles[0]);
    temp_file1 = working_dir + os.sep + "temp_dfj304jfd.txt";

    #Concatenate input files
    cat_command = ['cat'];
    cat_command.extend(inputfiles);
    with open(temp_file1, 'w') as fout:
        sp.check_call(cat_command, stdout=fout);

    #Sort file to be merged
    temp_file2 = working_dir + os.sep + "temp_fje094j3.txt";
    with open(temp_file2, 'w') as fout:
        sp.check_call(['sortBed','-i',temp_file1], stdout=fout);

    #Merge file
    if(output_file.find(os.sep) == -1):
        output_file = working_dir + os.sep + output_file;

    with open(output_file, 'w') as fout:
        sp.check_call(['bedtools','merge','-i',temp_file2], stdout=fout);

    #Clean up temporary files
    os.remove(temp_file1);
    os.remove(temp_file2);

    return output_file;

def bed_3cols(bedfile_name, outputfile_name=''):
    """
    Reads in <bedfile_name> and outputs a file with only the first
        three columns of each line
    """
    if(len(outputfile_name) == 0):
        outputfile_name = bedfile_name;

    with open(bedfile_name, 'r') as fin:
        contents = fin.read();

    csplit = contents.splitlines();


    with open(outputfile_name,'w') as fout:
        for i,line in enumerate(csplit):
            fout.write('\t'.join(line.split('\t')[0:3]) + '\n');

def submit_queue(run_directory, command):
    with open("run_template.sh", "r") as fin:
        template = fin.read();

    template = template.replace("<#LOCATION>", run_directory);
    template = template.replace("<#COMMAND>", command);

    outfile = run_directory + os.sep + "run.sh";

    with open(outfile, "w") as fout:
        fout.write(template);

    sp.check_call(["qsub", outfile]);

def download_data(antibody1, time1, antibody2, time2):

    exp_file = './Data/GSE36104-GPL15103_series_matrix.txt';
    exp_dict = read_series_matrix(exp_file);

    fastq_file1 = download_sra_for_sample(exp_dict, antibody1, time1);
    peaks_file1 = download_peaks_for_sample(exp_dict, antibody1, time1);

    fastq_file2 = download_sra_for_sample(exp_dict, antibody2, time2);
    peaks_file2 = download_peaks_for_sample(exp_dict, antibody2, time2);

    return (fastq_file1, peaks_file1, fastq_file2, peaks_file2);


def run_pipeline(antibody1, time1, antibody2, time2):

    working_dir = os.getcwd();

    #location of the experiment metadata file
    exp_file = meta_data.EXP_FILE_MATRIX;

    exp_dict = meta_data.read_series_matrix(exp_file);

    fastq_file1 = download_sra_for_sample(exp_dict, antibody1, time1);
    peaks_file1 = download_peaks_for_sample(exp_dict, antibody1, time1);

    fastq_file2 = download_sra_for_sample(exp_dict, antibody2, time2);
    peaks_file2 = download_peaks_for_sample(exp_dict, antibody2, time2);

    import process_raw;

    #Then can run ChipSeq pipeline on files
    bam1 = process_raw.ChIPSeq_Pipeline(fastq_file1);
    bed1 = bam_to_bed(bam1);

    bam2 = process_raw.ChIPSeq_Pipeline(fastq_file2);
    bed2 = bam_to_bed(bam2);

    bed_3cols(peaks_file1);
    bed_3cols(peaks_file2);
    peaks_filename = antibody1 + '_' + str(time1) + '_' + antibody2 + '_' + str(time2) + '.txt';
    peaks_file = bed_merge(peaks_filename, peaks_file1, peaks_file2);


    #MANorm

    from manorm import MANorm;

    manorm_analysis = MANorm();
    manorm_analysis.set_rawfiles(bed1, bed2);
    manorm_analysis.set_peakfiles(peaks_file1, peaks_file2);

    manorm_outputdir = working_dir + os.sep + "MANorm_Output";
    manorm_workdir = manorm_outputdir + "work_dir";
    manorm_analysis.set_workdir(manorm_workdir);
    manorm_analysis.set_outputdir(manorm_outputdir);

    manorm_analysis.run();

