from __future__ import division
import re
import editdistance
import time
import sys
import argparse

similarity_cutoff=0.1

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--Parsed_IgBlast_output', type=str)
parser.add_argument('-o', '--output_file', type=str)
args = parser.parse_args()
infile = args.Parsed_IgBlast_output
outfile = args.output_file
out=open(outfile,'w')


def collect_reads(sample_file):
    comb=[]
    compare={}
    for line in open(sample_file):
        a=line.strip('\n').split('\t')
        if len(a)>3:
            name,J,V,CDR3,isotype=a[0],a[4],a[2],a[5],a[8]
            if 0<len(CDR3)<200 and isotype!='N':
                 comb.append((name,V,J,CDR3))
                 if not compare.get(V+'_'+J):
                      compare[V+'_'+J]=[]
                 compare[V+'_'+J].append(line)
    return comb,compare

def cluster_reads(reads,compare):
    allmatches={}
    lineage_id=0

    for items in reads:
        matches=[]
        name,V,J,CDR3=items[0],items[1],items[2],items[3]
        cutoff=len(CDR3)*similarity_cutoff
        if not allmatches.get(items):
            lineage_id+=1
            allmatches[items]=lineage_id
            matches.append(items)
            new=True
            while new:
                new=False
                new_matches=[]
                counter=0
                for line in compare[V+'_'+J]:
                    a=line.strip('\n').split('\t')
                    name_c,J_c,V_c,CDR3_c=a[0],a[4],a[2],a[5]
                    items_c=(name_c,V_c,J_c,CDR3_c)
                    if not allmatches.get(items_c):
                        for entries in matches:
                            distance=editdistance.eval(CDR3_c,CDR3)
                            if distance<=cutoff:
#                                print('added',lineage_id,name_c)
                                allmatches[items_c]=lineage_id
                                new_matches.append(items_c)
                                new=True
                                break

                for match in new_matches:
                    matches.append(match)
    return allmatches

def write(compare,allmatches):
    entry_list=[]
    for combination in compare:
        entries=compare[combination]
        for line in entries:
            a=line.strip('\n').split('\t')
            name_c,J_c,V_c,CDR3_c=a[0],a[4],a[2],a[5]
            items_c=(name_c,V_c,J_c,CDR3_c)
            entry_list.append((allmatches[items_c],line.strip('\n')))
    for entry in sorted(entry_list, key= lambda x: x[0]):
        out.write(str(entry[0])+'\t'+entry[1]+'\n')


def main():
    print('collecting reads')
    reads,compare=collect_reads(infile)
    print('clustering reads')
    allmatches=cluster_reads(reads,compare)
    print('writing reads')
    write(compare,allmatches)


main()
