from __future__ import print_function, division;

import os;
from collections import namedtuple
import numpy as np;
import meta_data;
import re;
import sys;

nt = namedtuple("DiffAnalysis",["names", "scores", "pvalues", "sortkeys"]);

def num_from_chr(chr_string):
    root = chr_string[3:];
    if(root == "X"):
        return 20;
    if(root == "Y"):
        return 21;
    else:
        try:
            return int(root);
        except:
            raise ValueError("Unrecognized chr identifier " + chr_string);

def package_manorm(exp):
    """MaNorm data to common format"""
    files = meta_data.manorm_output(exp);

    for i,filename in enumerate(files):

        names = list();
        scores = list();
        pvalues = list();
        sortkeys = list();

        with open(filename, 'r') as fin:
            fin.readline();
            for line in fin:
                vals = line.strip().split('\t');
                start = int(float(vals[1]));
                stop = int(float(vals[2]));
                name = vals[0]+":" + str(start) + "-" + str(stop);
                score = float(vals[6]);
                pvalue = float(vals[8]);
                try:
                    sortkey = [num_from_chr(vals[0]), start, stop];
                except ValueError:
                    continue;  #Ignore things like ChrM
                
                names.append(name);
                scores.append(score);
                pvalues.append(pvalue);
                sortkeys.append(sortkey);

        scores = np.array(scores);
        pvalues = np.array(pvalues);
        pvalues = 10**(pvalues * -1);
        sortkeys = np.array(sortkeys);

        xx = nt(names, scores, pvalues, sortkeys);

        if(i == 0):
            unique_0 = xx;
        elif(i == 1):
            common_0 = xx;
        elif(i == 2):
            unique_X = xx;
        elif(i == 3):
            common_X = xx;
        else:
            raise Exception("Too many manorm files!!");

    out = (unique_0, common_0, unique_X, common_X);

    out = sort_package(out);

    return out;

def package_deseq_i(exp):
    """Gets deseq to an intermediate format
    Get to common format by aligning with MaNorm data
    """

    files = meta_data.deseq_output(exp);

    for i, filename in enumerate(files):

        names = list();
        scores = list();
        pvalues = list();
        sortkeys = list();

        with open(filename, 'r') as fin:
            fin.readline();
            for line in fin:
                vals = line.strip().split(',');
                start = int( re.search(":\d+-", vals[0]).group()[1:-1]);
                stop  = int( re.search("-\d+", vals[0]).group()[1:]);
                chrom = re.search("chr[^:]+", vals[0]).group();

                name = vals[0].replace('"', '');
                score = float(vals[2]);
                pvalue = float(vals[5]);
                try:
                    sortkey = [num_from_chr(chrom), start, stop];
                except ValueError:
                    continue;  #Ignore chrM for example

                names.append(name);
                scores.append(score);
                pvalues.append(pvalue);
                sortkeys.append(sortkey);

        scores = np.array(scores);
        pvalues = np.array(pvalues);
        sortkeys = np.array(sortkeys);

        xx = nt(names, scores, pvalues, sortkeys);

        if(i == 0):
            t0 = xx;
        elif(i == 1):
            tX = xx;
        else:
            raise Exception("Too many deseq outputs!!");

    return (t0, tX);

def package_foldchange_i(exp):
    """Same as the similar deseq method only for foldchange
    """
    files = meta_data.foldchange_output(exp);

    for i, filename in enumerate(files):

        names = list();
        scores = list();
        sortkeys = list();
        pvalues = list();

        with open(filename, 'r') as fin:
            fin.readline();
            for line in fin:
                vals = line.strip().split('\t');
                start = int(vals[1]);
                stop  = int(vals[2]);
                chrom = vals[0];
                pvalue = float(vals[4]);

                name = vals[0]+":"+str(start)+"-"+str(stop);
                score = float(vals[3]);
                try:
                    sortkey = [num_from_chr(chrom), start, stop];
                except ValueError:
                    continue;  #Ignore chrM for example

                names.append(name);
                scores.append(score);
                sortkeys.append(sortkey);
                pvalues.append(pvalue);

        scores = np.array(scores);
        sortkeys = np.array(sortkeys);
        pvalues = np.array(pvalues);

        #Generate pseudo p_values
        scores[scores == -1] = 1e3;
        scores[scores == 0] = 1e-3;

        #score_mag = np.max(np.vstack((scores, scores**-1)), axis=0);
        #pvalues = (score_mag.shape[0] - np.argsort(score_mag)) / score_mag.shape[0];
        
        xx = nt(names, scores, pvalues, sortkeys);

        if(i == 0):
            t0 = xx;
        elif(i == 1):
            tX = xx;
        else:
            raise Exception("Too many foldchange outputs!!");

    return (t0, tX);


def align_to_manorm(intermediate_out, manorm_out):
    """aligns intermediate deseq and foldchange outputs to manorm data
    matches based on sortkey sortings
    """

    #Split intermediate t0 into t0_unique and t0_common
    t0 = intermediate_out[0];
    ma_t0_unique = manorm_out[0];
    ma_t0_common = manorm_out[1];
    un_i = [i for i,name in enumerate(t0.names) if name in ma_t0_unique.names];
    comm_i = [i for i,name in enumerate(t0.names) if name in ma_t0_common.names];

    unique_t0_names = [t0.names[i] for i in un_i];
    unique_t0 = nt(unique_t0_names, t0.scores[un_i], t0.pvalues[un_i], t0.sortkeys[un_i,:]);

    common_t0_names = [t0.names[i] for i in comm_i];
    common_t0 = nt(common_t0_names, t0.scores[comm_i], t0.pvalues[comm_i], t0.sortkeys[comm_i,:]);

    #Split intermediate tX into tX_unique and tX_common
    tX = intermediate_out[1];
    ma_tX_unique = manorm_out[2];
    ma_tX_common = manorm_out[3];

    un_i = [i for i,name in enumerate(tX.names) if name in ma_tX_unique.names];
    comm_i = [i for i,name in enumerate(tX.names) if name in ma_tX_common.names];

    unique_tX_names = [tX.names[i] for i in un_i];
    unique_tX = nt(unique_tX_names, tX.scores[un_i], tX.pvalues[un_i], tX.sortkeys[un_i,:]);

    common_tX_names = [tX.names[i] for i in comm_i];
    common_tX = nt(common_tX_names, tX.scores[comm_i], tX.pvalues[comm_i], tX.sortkeys[comm_i,:]);

    return (unique_t0, common_t0, unique_tX, common_tX);

def package_deseq(exp, manorm_out):

    intermediate = package_deseq_i(exp);
    out = align_to_manorm(intermediate, manorm_out);

    out = sort_package(out);

    return out;

def package_foldchange(exp, manorm_out):

    intermediate = package_foldchange_i(exp);
    out = align_to_manorm(intermediate, manorm_out);

    out = sort_package(out);

    return out;

def sort_package(package_data):
  
    out_list = list();

    for data in package_data:
        sortkeys = data.sortkeys;

        sort_records = np.core.records.fromrecords(sortkeys, names="chr,start,stop");

        ii = np.argsort(sort_records, order=["chr", "start", "stop"]);

        ordered_names = [data.names[i] for i in ii];
        ordered_scores = data.scores[ii];
        ordered_pvalues = data.pvalues[ii];
        ordered_sortkeys = data.sortkeys[ii,:];

        out_list.append(nt(ordered_names, ordered_scores, ordered_pvalues, ordered_sortkeys));


    return tuple(out_list);

def ma_label_peaks(ma_data, p_threshold=0.05, ratio_threshold = 1.0):
    """labels each peak as either 1 (upregulatied), -1 (downregulated), or 0 (no change)
    output vector is a concatenation of peaks in t0_unique, t0_common, and tX_unique
    """
    ratio_threshold = np.log2(ratio_threshold)
    all_names = ma_data[0].names + ma_data[1].names + ma_data[3].names;
    all_scores = np.hstack((ma_data[0].scores, ma_data[1].scores, ma_data[3].scores));
    all_pvals = np.hstack((ma_data[0].pvalues, ma_data[1].pvalues, ma_data[3].pvalues));

    all_scores = all_scores * -1;  #MANorm does peak1/peak2, we want peak2/peak1

    all_labels = np.zeros(all_scores.shape);

    upreg_i =   np.logical_and(all_pvals <= p_threshold, all_scores >= ratio_threshold);
    downreg_i = np.logical_and(all_pvals <= p_threshold, all_scores <= -1*ratio_threshold);

    all_labels[upreg_i] = 1;
    all_labels[downreg_i] = -1;

    return all_labels;


def de_label_peaks(de_data, p_threshold = 0.05):
    """labels each peak as either 1 (upregulatied), -1 (downregulated), or 0 (no change)
    output vector is a concatenation of peaks in t0_unique, t0_common, and tX_unique
    """
    all_scores = np.hstack((de_data[0].scores, de_data[1].scores, de_data[3].scores));
    all_pvals = np.hstack((de_data[0].pvalues, de_data[1].pvalues, de_data[3].pvalues));

    #Wow, deseq p-values basically show nothing is significant ever.
    #So...I'm just going to take the top 5% peaks like in fold-change

    p_i = np.argsort(all_pvals);
    thresh_p = all_pvals[p_i[np.ceil(all_pvals.shape[0] * p_threshold)]];  

    all_labels = np.zeros(all_scores.shape);

    upreg_i =   np.logical_and(all_pvals <= thresh_p, all_scores > 0);
    downreg_i = np.logical_and(all_pvals <= thresh_p, all_scores < 0);

    all_labels[upreg_i] = 1;
    all_labels[downreg_i] = -1;

    return all_labels;

def fc_label_peaks(fc_data, p_threshold = 0.05, ratio_threshold = 1.0):
    """labels each peak as either 1 (upregulatied), -1 (downregulated), or 0 (no change)
    output vector is a concatenation of peaks in t0_unique, t0_common, and tX_unique
    """
    all_names = fc_data[0].names + fc_data[1].names + fc_data[2].names;
    all_scores = np.hstack((fc_data[0].scores, fc_data[1].scores, fc_data[3].scores));
    all_pvals = np.hstack((fc_data[0].pvalues, fc_data[1].pvalues, fc_data[3].pvalues));
    
    #Start new region
    #all_scorex = np.abs(np.log2(all_scores));
    #yy = np.argsort(all_scorex)[::-1];
    #good_p = yy[0:np.floor(len(all_scorex)*threshold)];
    #all_pvals[:] =1;
    #all_pvals[good_p] = 0;
    #End modified region

    all_labels = np.zeros(all_scores.shape);

    upreg_i =   np.logical_and(all_pvals <= p_threshold, all_scores >= ratio_threshold);
    downreg_i = np.logical_and(all_pvals <= p_threshold, all_scores <= 1/ratio_threshold);

    all_labels[upreg_i] = 1;
    all_labels[downreg_i] = -1;

    return all_labels;

def overlaps_vs_p(exp):
    ma_data = package_manorm(exp);
    fc_data = package_foldchange(exp, ma_data);

    p_values = 10**(np.array([0.0, -1.0, -2, -3, -4, -5, -6, -7, -8, -9, -10]));
    number_peaks_ma = np.zeros(p_values.shape);
    number_peaks_fc = np.zeros(p_values.shape);
    number_peaks_overlap = np.zeros(p_values.shape);


    for i in range(len(p_values)):
        ma_labels = ma_label_peaks(ma_data, p_values[i]);
        fc_labels = fc_label_peaks(fc_data, p_values[i]);

        number_peaks_ma[i] = np.count_nonzero(ma_labels != 0);
        number_peaks_fc[i] = np.count_nonzero(fc_labels != 0);

        number_peaks_overlap[i] = np.count_nonzero(np.logical_and(
            ma_labels != 0,
            ma_labels == fc_labels));

    return p_values, number_peaks_ma, number_peaks_fc, number_peaks_overlap;

def post_process_exp(exp):

    ma_data = package_manorm(exp);
    de_data = package_deseq(exp, ma_data);
    fc_data = package_foldchange(exp, ma_data);

    #Analyze results here

    #Return vector denoting each peak as UP, DOWN, or No Change (1, -1 o r 0)
    ma_labels = ma_label_peaks(ma_data, p_threshold = 1e-10); 
    de_labels = de_label_peaks(de_data, p_threshold = 1e-10);
    fc_labels = fc_label_peaks(fc_data, p_threshold = 1e-10);

    #Compute overlaps on that (vector of 7 denoting # of equal elements)
    overlaps = [0]*7;

    overlaps[0] = np.count_nonzero(ma_labels != 0);
    overlaps[1] = np.count_nonzero(de_labels != 0);
    overlaps[2] = np.count_nonzero(fc_labels != 0);

    overlaps[3] = np.count_nonzero(np.logical_and(
        ma_labels != 0,
        ma_labels == de_labels));

    overlaps[4] = np.count_nonzero(np.logical_and(
        ma_labels != 0,
        ma_labels == fc_labels));

    overlaps[5] = np.count_nonzero(np.logical_and(
        de_labels != 0,
        de_labels == fc_labels));

    overlaps[6] = np.count_nonzero(np.logical_and(
        ma_labels != 0,
        np.logical_and(ma_labels == de_labels, ma_labels==fc_labels)));

    overlaps = np.array(overlaps);

    #Summarize up-down counts
    ma_summary = [np.count_nonzero(ma_labels == 1),
                  np.count_nonzero(ma_labels == -1),
                  np.count_nonzero(ma_labels == 0)];

    de_summary = [np.count_nonzero(de_labels == 1),
                  np.count_nonzero(de_labels == -1),
                  np.count_nonzero(de_labels == 0)];

    fc_summary = [np.count_nonzero(fc_labels == 1),
                  np.count_nonzero(fc_labels == -1),
                  np.count_nonzero(fc_labels == 0)];

    summary = [ma_summary, de_summary, fc_summary];
    summary = np.array(summary);


    #Find potential interesting regions?
    interesting_name = list();
    interesting_score = list();

    for x,y in zip(ma_data, de_data):
        xi = np.argsort(x.scores);
        yi = np.argsort(y.scores);

        dd = np.abs(xi - yi);
        di = np.argsort(dd)[::-1]; #Sort differences descending

        if(len(di) > 10):
            di = di[0:10];

        interesting_name.extend([x.names[i] for i in di]);
        interesting_score.extend([x.scores[i] for i in di]);

    for x,y in zip(ma_data, fc_data):
        xi = np.argsort(x.scores);
        yi = np.argsort(y.scores);

        dd = np.abs(xi - yi);
        di = np.argsort(dd)[::-1]; #Sort differences descending

        if(len(di) > 10):
            di = di[0:10];

        interesting_name.extend([x.names[i] for i in di]);
        interesting_score.extend([x.scores[i] for i in di]);

    for x,y in zip(de_data, fc_data):
        xi = np.argsort(x.scores);
        yi = np.argsort(y.scores);

        dd = np.abs(xi - yi);
        di = np.argsort(dd)[::-1]; #Sort differences descending

        if(len(di) > 10):
            di = di[0:10];

        interesting_name.extend([x.names[i] for i in di]);
        interesting_score.extend([x.scores[i] for i in di]);



    #Take top 10 "interesting" sites

    i_score = np.array(interesting_score);

    ii = np.argsort(i_score)[::-1];

    interesting_name = [interesting_name[i] for i in ii[0:10]];
    interesting_score = [interesting_score[i] for i in ii[0:10]];

    #Save results
    out_directory = meta_data.directory_results(exp);

    if(not os.path.isdir(out_directory)):
        os.mkdir(out_directory);

    np.save(meta_data.overlaps_file(exp),overlaps);
    np.save(meta_data.summary_file(exp), summary);

    with open(meta_data.interesting_file(exp),'w') as fout:
        for name, score in zip(interesting_name, interesting_score):
            chrom = re.search(".*:", name).group()[:-1];
            start = re.search(":\d+", name).group()[1:];
            stop = re.search("-\d+", name).group()[1:];
            fout.write("\t".join([chrom, start, stop, str(score)]) + "\n");


def aggregate_overlaps(diff):
    overlaps = np.zeros(7);

    for pair in diff:
        o_file = meta_data.overlaps_file(pair[1]);

        o = np.load(o_file);
        overlaps = overlaps + o;

    return overlaps; 

def aggregate_summaries(diff):
    summaries = np.zeros((3,3));

    for pair in diff:
        s_file = meta_data.summary_file(pair[1]);

        s = np.load(s_file);
        summaries = summaries + s;

    return summaries;


if __name__ == "__main__":
    directory = sys.argv[1];
    exp = meta_data.exp_from_directory(directory);
    post_process_exp(exp);
