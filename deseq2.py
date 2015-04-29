from __future__ import print_function, division;

import os;
import sys;
import meta_data;
import subprocess as sp;
import file_operations;
import re;
from pybedtools import BedTool;

LOCATION = os.path.dirname(os.path.abspath(__file__));
R_SOURCE_FILE = "deseq2_script.r";
INPUT_FILE_NAME = "rpy_input.csv";

def generate_inputs(peaks, cond1, cond2, directory):
    peaks_bed = BedTool(peaks);
    cond1_bed = BedTool(cond1);
    cond2_bed = BedTool(cond2);

    over1 = peaks_bed.intersect(cond1_bed, c=True);
    over2 = peaks_bed.intersect(cond2_bed, c=True);

    lines = list();  #Horrible way of doing it
    #But for some reason, pybedtools crashes if I try to iterate through
    #a bedtool while writing a file

    #Write a csv describing the peaks for deseq2
    for a,b in zip(over1, over2):
        name = a.chrom + ":" + str(a.start) + "-" + str(a.stop);
        line = ",".join([name, str(a.count), str(b.count)]) + '\n';
        lines.append(line);

    out_file = directory + os.sep + INPUT_FILE_NAME;
    with open(out_file,"w") as fout:
        fout.write("untreated" + "," + "treated" + "\n");
        for line in lines:
            fout.write(line);

def setup_r_script(directory):
    with open(LOCATION + os.sep + R_SOURCE_FILE) as fin:
        source_file = fin.read();

    source_file = source_file.replace("<#COUNT_CSV>", directory + os.sep + INPUT_FILE_NAME);

    with open(directory + os.sep + R_SOURCE_FILE, 'w') as fout:
        fout.write(source_file);

def deseq_analysis(exp1, exp2):
    #Gather files
    rawdir1 = meta_data.directory_raw_reads(exp1);
    rawdir2 = meta_data.directory_raw_reads(exp2);

    files = os.listdir(rawdir1);
    for filename in files:
        if(filename.endswith(".bed")):
            counts1 = rawdir1 + os.sep + filename;

    files = os.listdir(rawdir2);
    for filename in files:
        if(filename.endswith(".bed")):
            counts2 = rawdir2 + os.sep + filename;

    peaks = meta_data.peaks_file(exp1);

    out_directory = meta_data.directory_deseq(exp2);

    if(not os.path.isdir(out_directory)):
        os.mkdir(out_directory);

    generate_inputs(peaks, counts1, counts2, out_directory);
    setup_r_script(out_directory);

    #Call R script
    sp.check_call(['Rscript', out_directory + os.sep + R_SOURCE_FILE]);


if __name__ == "__main__":
    antibody = sys.argv[1];
    time = sys.argv[2];

    for exp in meta_data.exp_list:
        etime = int(re.search("\d+",exp["time"]).group());
        if(exp["chip antibody"] == antibody):
            if(etime == 0):
                exp1 = exp;
            if(etime == time):
                exp2 = exp;

    deseq_analysis(exp1, exp2);



