from __future__ import division, print_function

import os
import sys
import pybedtools
from pybedtools import BedTool;
import file_operations
import meta_data

def compute_fold_change(exp1, exp2):
    peaks_file = meta_data.peaks_file(exp1);

    raw_file1 = meta_data.raw_bed_file(exp1);
    raw_file2 = meta_data.raw_bed_file(exp2);

    peaks = BedTool(peaks_file);
    raw1 = BedTool(raw_file1);
    raw2 = BedTool(raw_file2);

    coverage_1 = peaks.intersect(raw1, c=True);
    coverage_2 = peaks.intersect(raw2, c=True);

    output_bed = list();
    #Bad way to do this, but I'm having trouble writing
    #to a file while iterating over BedTools

    for i1, i2 in zip(coverage_1, coverage_2):
        if(i1.count == 0):
            cstring = "-1";
        else:
            cstring = str.format("{:f}",i2.count / i1.count);
        line = i1.chrom + "\t" + str(i1.start) + "\t" + str(i1.end) + "\t" + cstring + "\n";
        output_bed.append(line);

    out_directory = meta_data.directory_foldchange(exp2);

    if(not os.path.isdir(out_directory)):
        os.mkdir(out_directory);

    out_file = out_directory + os.sep + "foldchange.bed";

    with open(out_file, 'w') as fout:
        for line in output_bed:
            fout.write(line);


if __name__ == "__main__":
    dir1 = sys.argv[1];
    dir2 = sys.argv[2];
    
    exp1 = meta_data.exp_from_directory(dir1);
    exp2 = meta_data.exp_from_directory(dir2);

    compute_fold_change(exp1, exp2);
