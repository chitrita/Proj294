"""
Wrapper file to call macs2 in python
"""

import meta_data;
import os;
import sys;
import subprocess as sp;

def macs2_callpeaks(exp):
    """Calls peaks on the experiments raw file.
    Puts the result in the experiments Peaks folder.
    """

    input_file = meta_data.raw_bed_file(exp);

    name = "macs2_" + meta_data.exp_antibody(exp) + "_" + str(meta_data.exp_time(exp));

    outdir = meta_data.directory_peaks(exp);

    g = "1.87e9";

    sp.check_call(["macs2", "callpeak", "-t", input_file, "-f", "BED", "-n", name, "--outdir", outdir, "-g", g]);

if __name__ == "__main__":
    directory = sys.argv[1];

    exp = meta_data.exp_from_directory(directory);

    macs2_callpeaks(exp);
