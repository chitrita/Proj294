from pybedtools import BedTool;
import numpy as np;


def generate_curve(bedfile, chromosome, region_start, region_stop):
    bedtool = BedTool(bedfile);
    region_of_interest = BedTool(chromosome + ' ' + str(region_start) + ' ' + str(region_stop), from_string=True);

    plot_region = region_of_interest.intersect(bedtool);

    domain = np.arange(region_start, region_stop+1);

    values = np.zeros(domain.shape);

    for interval in plot_region:
        if(interval.start < region_start):
            start = region_start;
        else:
            start = interval.start;

        if(interval.end > region_stop):
            finish = region_stop;
        else:
            finish = interval.end;

        start_i = start - domain[0];
        finish_i = finish - domain[-1];

        values[start_i:finish_i] = values[start_i:finish_i] + 1;

    return (domain, values);



