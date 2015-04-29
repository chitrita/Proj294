from __future__ import print_function, division;
import os;
import re;

ROOTDIR = "/data/yosef/users/david.detomaso/Proj294";

CONTROL_FILE_MATRIX = "/home/eecs/david.detomaso/Proj294/Data/GSE36104-GPL11002_series_matrix.txt";

EXP_FILE_MATRIX = "/home/eecs/david.detomaso/Proj294/Data/GSE36104-GPL15103_series_matrix.txt";

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
        for j in range(1,10):
            if(key == ('Sample_supplementary_file_' + str(j))):
                exp_dict[key] = split_line[1:];

    return exp_dict;


def select_exp(exp_dict, antibody, time, all=False):

    match_antibody = {i for i, x in enumerate(exp_dict['chip antibody']) if antibody.lower() == x.lower()}

    corrected_time = str(time) + " min";

    match_time = {i for i, x in enumerate(exp_dict['time']) if corrected_time.lower() == x.lower()}

    match_sample = match_antibody.intersection(match_time);
    if(len(match_sample) == 0):
        print("Error, no matching sample");
        return;
    if(len(match_sample) > 1):
        if(not all):
            print("Error, multiple matching samples");
            return;
        else:
            all_exps = list();
            for sample in match_sample:
                exp_info = dict();
                for key in exp_dict.keys():
                    exp_info[key] = exp_dict[key][sample];
                all_exps.append(exp_info);

            return all_exps;
    else:  #length equals one

        match_sample = match_sample.pop();

        exp_info = dict();
        for key in exp_dict.keys():
            exp_info[key] = exp_dict[key][match_sample];

        return exp_info;



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

def exp_from_directory(directory):
    for exp in exp_list:
        if(directory_from_exp(exp) == directory):
            return exp;
    raise Exception("No Matching Experiment!");

def directory_raw_reads(exp):
    return directory_from_exp(exp) + os.sep + "RAW";

def directory_peaks(exp):
    return directory_from_exp(exp) + os.sep + "Peaks";

def directory_manorm(exp):
    return directory_from_exp(exp) + os.sep + "MANorm_MACS";

def directory_deseq(exp):
    return directory_from_exp(exp) + os.sep + "DeSeq";

def directory_foldchange(exp):
    return directory_from_exp(exp) + os.sep + "FoldChange";

def peaks_file(exp):
    """Returns the location of the peaks file for a given experiment"""
    directory = directory_peaks(exp);
    files = os.listdir(directory);
    #peaks_re = re.compile("GSM.*peaks.*\.txt");  #Downloaded Peaks File
    peaks_re = re.compile(".narrowPeak");
    for filename in files:
        if(peaks_re.search(filename) is not None):
            return directory + os.sep + filename;

def raw_file(exp):
    """Returns the location of the raw file for a given experiment"""
    directory = directory_raw_reads(exp);
    files = os.listdir(directory);
    for filename in files:
        if(filename == "aligned_sorted_marked_duplicates_" + str(exp_time(exp)) + ".bam"):
            return directory + os.sep + filename;
    raise Exception("Raw File not Found");

def raw_bed_file(exp):
    directory = directory_raw_reads(exp);
    for filename in os.listdir(directory):
        if(filename.endswith(".bed")):
            return directory + os.sep + filename;
    raise OSError("Raw Bed file not found");

def exp_time(exp):
    return(int(re.search("\d+", exp["time"]).group()));

def exp_antibody(exp):
    return exp["chip antibody"];

def define_diff_pairs():
    """Returns a list of 2-tuples each with one exp object.
    Each pair is used in a differential analysis.  All are Time X vs Time 0.
    Results stored in Time X's folder.
    """
    diff_pairs = list();
    int_reg = re.compile('\d+');
    for exp in exp_list:
        time = int(int_reg.search(exp['time']).group());
        antibody = exp['chip antibody'];
        if(time > 0):
            control_exp = select_exp(exp_dict, antibody, 0);
            diff_pairs.append((control_exp, exp));

    return diff_pairs;

