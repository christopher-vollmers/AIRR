import os
import sys
import re
import editdistance
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--IgBlast_file', type=str)
parser.add_argument('-t', '--isotype_fasta', type=str)
parser.add_argument('-f', '--isoform_fasta', type=str)
parser.add_argument('-o', '--output_file', type=str)
parser.add_argument('-s', '--test', action='store_true',default=False)



args = parser.parse_args()
IgBlast_file = args.IgBlast_file
isotype_fasta = args.isotype_fasta
isoform_fasta = args.isoform_fasta
output_file = args.output_file
test=args.test

handle = open(IgBlast_file)
out=open(output_file,'w')
header = handle.readline()


def determine_isoform(C_seq,match,IsoDict):
    typeDict={}
    typeDict['S']=0
    typeDict['M']=0
    for i in range(0,len(C_seq),3):
        kmer=C_seq[i:i+12]
        if kmer in IsoDict:
            matches=IsoDict[kmer]
            for isoform in matches:
                if match in isoform:
                     type1=isoform.split('_')[1]
                     typeDict[type1]+=1
    if typeDict['S']>typeDict['M']:
         isoform='S'

    elif typeDict['S']<typeDict['M']:
         isoform='M'
    else:
         isoform='-'
    return isoform,typeDict['S'],typeDict['M']

def find_constant_region(C_seq,Isotypes):

    match='_'
    min1=10000000

    for Isotype in Isotypes:

        C_limit=min(len(C_seq),len(Isotype[1]),300)
        if C_limit>8:
            distance1=editdistance.eval(C_seq[:C_limit],Isotype[1][:C_limit])
            if distance1/(C_limit+1)<0.15:
                if distance1<min1:
                    min1=distance1
                    match=Isotype[0]
    return match,min1

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readList = []
    tempSeqs, headers, sequences = [], [], []
    for line in open(inFile):

        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            headers.append(line.split()[0][1:])
        # covers the case where the file ends while reading sequences
        if line.startswith('>'):
            sequences.append(''.join(tempSeqs).upper())
            tempSeqs = []
        else:
            tempSeqs.append(line)
    sequences.append(''.join(tempSeqs).upper())
    sequences = sequences[1:]
    for i in range(len(headers)):
        readList.append((headers[i],sequences[i]))
    return readList

def read_constant_region_annotations(isoform_fasta):
    Isoforms=read_fasta(isoform_fasta)
    IsoDict={}
    for isoform,sequence in Isoforms:
        if isoform not in IsoDict:
            for i in range(0,120,1):
                if sequence[i:i+12] not in IsoDict:
                    IsoDict[sequence[i:i+12]]=[]
                IsoDict[sequence[i:i+12]].append(isoform)
    return IsoDict

def process_reads(handle,Isotypes,IsoDict,test):
    header_list = header.strip().split("\t")
    for line in handle:
        la = line.strip().split("\t")
        if len(la) == 88:
            V_seg,D_seg,J_seg = la[7].split(",")[0],la[8].split(",")[0],la[9].split(",")[0]
            if len(la[7].split(","))>1:
                Amb='Ambiguous:'+la[7]
            else:
                Amb='N/A'
            germline_mismatch_pos = []
            sequence_V,germline_V,CDR3,seq_type,ID,seq = la[20],la[22],la[42],la[2],la[0],la[1]
            s_index,g_index,C_seq = int(la[60]),int(la[62]),seq[int(la[69])-1:-20]

            match,distance=find_constant_region(C_seq,Isotypes)
            if test:
                isoform,S,M=determine_isoform(C_seq,match,IsoDict)
#                print(match,distance,isoform,S,M)
            else:
                isoform='N/A'

            for s,g in zip(sequence_V,germline_V):
                if s != '-' and g != '-' and s != g:
                    germline_mismatch_pos.append(str(g_index))
                if s != "-":
                    s_index+=1
                if g != '-':
                    g_index+=1

            string='%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ID.replace('reversed|',''),seq_type,V_seg,D_seg,J_seg,CDR3,Amb,isoform,match,seq,",".join(germline_mismatch_pos))
            out.write(string)

def main():
    print('reading isotype sequences')
    Isotypes=read_fasta(isotype_fasta)
    if test:
        print('reading isoform sequences')
        IsoDict=read_constant_region_annotations(isoform_fasta)
    else:
        IsoDict={}
    print('processing reads')
    process_reads(handle,Isotypes,IsoDict,test)

main()
