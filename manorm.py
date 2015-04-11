# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 22:38:12 2015

@author: daved_000
"""
from __future__ import print_function, division;
import os;
import subprocess as sp;

matlab_source = """

addpath('CodeLib');

%%%%%% start of parameter region %%%%%%
nchr = 21;              % number of chromosomes; 24 for human, 21 for mouse

% the bed file name of peak coordinate of sample 1
peakfile_1 = '$peakfile_1$';
% the bed file name of read coordinate of sample 1: Yale ENCODE MYC ChIP-seq data replicate 1 in HelaS3 cells
rawdatafile_1 = '$rawdatafile_1$';
% shift size of reads for sample 1; which equal to the half of the DNA fragment after size selection
tag_shift_1 = 100;      
% length of reads for sample 1
tag_length_1 = 36;        

% the bed file name of peak coordinate of sample 2
peakfile_2 = '$peakfile_2$;
% the bed file name of read coordinate of sample 2: Yale ENCODE MYC ChIP-seq data replicate 1 in K562 cells
rawdatafile_2 = '$rawdatafile_2$';
% shift size of reads for sample 1; which equal to the half of the DNA fragment after size selection
tag_shift_2 = 100;      
% length of reads for sample 1
tag_length_2 = 36;   

% directory of work folder to store temperary files; should be made before analysis
workdir = '$workdir$';
% directory of output folder to store output files; should be made before analysis
outputdir = '$outputdir$'; 
% name of comparison
MApair = 'Experiment Name';

extention = 500;            % 2*extention = the window size to calculate read density; we recommend extention=1000 for histone modification, extention=500 for transcription factor
summit2summit_dist = 250;      % the selected common peaks much have summit-to-summit distance smaller than this value to as matched common peaks; suggest to be half of extention

Mcut_unbiased = 1;          % Curoff of |M| value to define unbiased peaks between 2 samples; here M value is the log2 fold change at peak region between two samples, value 1 means fold change 2
                            % Non-cell type-specific binding regions: all peaks associated with |M|<Mcut_unbiaseda.
                            
Mcut_biased = 1;            % Cutoff of M value cutoff to define biased peaks between 2 samples: 
Pcut_biased = 1e-2;         % Cutoff of P value (of differential binding) cutoff to define biased peaks between 2 samples;
                            % Sample 1-specific binding regions: all peaks associated with M>=Mcut_biased and P<Pcut_biased; 
                            % Sample 2-specific binding regions: all peaks associated with M<=-Mcut_biased and P<Pcut_biased;
%%%%%% end of parameter region %%%%%%

$OTHERPARAMS$

% read peaks;
disp(' ');
disp('Step1: read peaks');
peak_name_1 = step1_read_MACS_peaks(peakfile_1, nchr, workdir);
peak_name_2 = step1_read_MACS_peaks(peakfile_2, nchr, workdir);

% classify common or unique peaks
disp(' ');
disp('Step2: classify common or unique peaks');
Nsamp = 100;                 % number of random permutation to calculate fold enrichment of overlap 
step2_classify_2peak_sets_by_overlap(peak_name_1, peak_name_2, workdir, Nsamp);

% split raw data by chromosome
disp(' ');
disp('Step3: read and split raw data by chromosome');
[pathstr1, marker_name_1, ext1, versn1] = fileparts(rawdatafile_1);
step3_split_sequencing_data_2chr(rawdatafile_1, nchr, workdir);

[pathstr2, marker_name_2, ext2, versn2] = fileparts(rawdatafile_2);
step3_split_sequencing_data_2chr(rawdatafile_2, nchr, workdir);

% calculate read density in peaks
disp(' ');
disp('Step4: calculate read density in peaks');
step4_calculate_peak_read_density(peak_name_1, peak_name_2, marker_name_1, marker_name_2, extention, workdir, nchr, tag_shift_1, tag_shift_2, tag_length_1, tag_length_2);

% calculate normalized M-A values 
disp(' ');
disp('Final Step: normalize and output the M-A values');
step5_rescale_peak_strength_by_common_peaks(peak_name_1, peak_name_2, marker_name_1, marker_name_2, summit2summit_dist, workdir, MApair, outputdir);

% output the normalized M-value to wig file (in order to upload to genome browsers)
disp(' ');
disp('Summarize data: output the M value of 2 peak sets to wig files');
step6_output_peak_Mvalue_2wig_file(peak_name_1, workdir, MApair, outputdir);
step6_output_peak_Mvalue_2wig_file(peak_name_2, workdir, MApair, outputdir);

% merge common peaks and define unbiased and biased peaks between two samples based on given M-value and P-value cutoff 
disp(' ');
disp('Summarize data: merge common peaks and extract differential and non-differential binding regions');
step7_define_biased_unbiased_peaks(peak_name_1, peak_name_2, workdir, MApair, outputdir, Mcut_unbiased, Mcut_biased, Pcut_biased)

% clean up
rmpath('CodeLib');
clear all;
""";

class MANorm:
    
    matlab_filename = 'matlab_script.m';    
    
    def __init__(self):
        self.matlab_source = matlab_source;
        self.peakfile1 = '';
        self.peakfile2 = '';
        self.rawfile1 = '';
        self.rawfile2 = '';
        self.workdir = '';
        self.outputdir = '';
        self.otherparams = [''];
        
    def run(self):
        
        if(len(self.peakfile1) == 0 or len(self.peakfile2) == 0):
            print("Error, Peak files not set");
            return;
        
        if(len(self.rawfile1) == 0 or len(self.rawfile2) == 0):
            print("Error, Raw files not set");
            return;
        
        if(len(self.workdir) == 0):
            print("Error, Work Directory not set");
            return;
            
        if(len(self.outputdir) == 0):
            print("Error, Output Directory not set");
            return;
            
        self.matlab_source = self.matlab_source.replace('$peakfile_1$',self.peakfile1);
        self.matlab_source = self.matlab_source.replace('$rawdatafile_1$',self.rawfile1);
        self.matlab_source = self.matlab_source.replace('$peakfile_2$',self.peakfile2);
        self.matlab_source = self.matlab_source.replace('$rawdatafile_2_1$',self.rawfile2);
        self.matlab_source = self.matlab_source.replace('$workdir$',self.workdir);
        self.matlab_source = self.matlab_source.replace('$outputdir$',self.outputdir);
        self.matlab_source = self.matlab_source.replace('$OTHERPARAMS$', ';\n'.join(self.otherparams));

        ff = open(self.workdir + os.sep + self.matlab_filename, 'w');
        ff.write(self.matlab_source);
        ff.close();
        
        sp.check_call(['matlab', '-r', self.matlab_filename]);
        
    def set_peakfiles(self, peakfile1, peakfile2):
        self.peakfile1 = peakfile1;
        self.peakfile2 = peakfile2;
    
    def set_rawfiles(self, rawfile1, rawfile2):
        self.rawfile1 = rawfile1;
        self.rawfile2 = rawfile2;
    
    def set_outputdir(self, outputdir):
        self.outputdir = outputdir;
        
    def set_workdir(self, workdir):
        self.workdir = workdir;
    
    