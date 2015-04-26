from __future__ import division, print_function;

import meta_data;
import file_operations;
import os;

def create_directories():
    for exp in meta_data.exp_list:
        directory = meta_data.directory_from_exp(exp);
        os.makedirs(directory);

def download_raw_reads():
    for i,exp in enumerate(meta_data.exp_list):
        print(i, exp["chip antibody"], exp["time"]);
        file_operations.download_sra_for_sample(exp);


