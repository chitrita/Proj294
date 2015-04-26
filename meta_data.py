from __future__ import print_function, division;
import os;
from file_operations import *;

ROOTDIR = "/data/yosef/users/david.detomaso/Proj294";


AA_LIST = ["Ctcf",
           "PU1",
           "Cebpb",
           "Junb",
           "Irf4",
           "Irf1",
           "Rela",
           "Egr2",
#           "Runx1",  has no 120
           "Maff",
           "E2f1",
           "Ahr",
           "Stat3",
#           "Relb",  these two have duplicated paired-end test
#           "Rel",
           "Stat1",
           "E2f4",
           "Stat2",
           "Nfkb1",
           "Egr1",
           "Irf2",
           "Hif1a"];

TIMES = [0, 30, 60, 120];

exp_dict = read_series_matrix(EXP_FILE_MATRIX);

exp_list = list();
for aa in AA_LIST:
    for time in TIMES:
        #print(aa, time);
        exp_list.append(select_exp(exp_dict, aa, time));

for aa in ["Runx1"]:
    for time in [0, 60, 120]:
        exp_list.append(select_exp(exp_dict, aa, time));

for aa in ["Relb", "Rel"]:
    for time in [0, 30, 60]:
        exp_list.append(select_exp(exp_dict, aa, time));

    time = 120;
    multi_exps = select_exp(exp_dict, aa, time, all=True);
    for exp in multi_exps:
        if(len(exp["chip antibody catalog number"]) > 0):  #the single ended exp has this property in the matrix
            exp_list.append(exp);

def directory_from_exp(exp):
    etime = exp["time"].replace(" ", "");
    return ROOTDIR + os.sep + exp["chip antibody"] + os.sep + etime;

def directory_raw_reads(exp):
    return directory_from_exp(exp) + os.sep + "RAW";

