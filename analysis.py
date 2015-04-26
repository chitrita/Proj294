from pybedtools import BedTool; 
from itertools import combinations;
from collections import namedtuple;

def compute_all_overlaps(*peak_calls):
    """
    Each peak_call in peak_calls should be a namedtuple with fields
    name, and file.

    file is a bed file defining peak regions

    Computes overlaps for every 2-pair, 3-group, 4-group, etc between input peak calls
    
    Returns a summary object
    """
    
    result = list();

    Group = namedtuple("Group", ["peak_calls", "overlap"]);

    for num in range(len(peak_calls)+1):
        for group in combinations(peak_calls, num):
            if(len(group)>1):
                count = compute_overlap(group);
                names = [x.name for x in group];
                result.append(Group(names, count));
            else #only one group member
                pbed = BedTool(group[0].file);
                count = pbed.count();
                result.append(Group(group[0].name, count));

    return result;



def compute_overlap(*peak_calls):
    """
    Computes the overlap between the supplied set of peak calls.
    Each peak_call in peak_calls is a namedtuple as above.
    """

    original = BedTools(peak_calls[0].file);

    for peak_call in peak_calls[1:]:
        original = original.intersect(BedTools(peak_call.file));

    count = original.count();

    return count;

