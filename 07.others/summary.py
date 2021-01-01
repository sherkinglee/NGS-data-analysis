#!/usr/bin/env python
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2020-12-31 18:49:06
LastEditors: Li Fajin
LastEditTime: 2021-01-01 16:42:11
Description: summary of processed reads from all steps
'''


import re
import os
import sys
from optparse import OptionParser


def create_parser_for_reads_summary():
	'''argument parser.'''
	usage="usage: python %prog "
	parser=OptionParser(usage=usage)
	parser.add_option("-c","--cutadapt",action="store",type="string",dest="cutadaptLog",
			help="Log file generated by cutadapt.")
	parser.add_option("-f","--filtering",action="store",type="string",dest="filteringLog",
			help="Log file generated by fastq_quality_filter.")
	parser.add_option("-r","--rRNA-contam",action="store",type="string",dest="rRNAcontamLog",
			help="Log file generated by rRNA contamination.")
	parser.add_option("-m","--star-mapping",action="store",type="string",dest="starMapping",
			help="Log file generated by star mapping.")
	parser.add_option("-s","--statistics",action="store",type="string",dest="statistics",
			help="Log file generated by STATISTICREADSONDNASCONTAM.")
	parser.add_option("--oc",action="store",type="string",dest="outCutadapt",help="output file for cutadapt log")
	parser.add_option("--of",action="store",type="string",dest="outFiltering",help="output file for filtering log")
	parser.add_option("--or",action="store",type="string",dest="outrRNAcontam",help="output file for rRNA contamination log")
	parser.add_option("--om",action="store",type="string",dest="outMapping",help="output file for star mapping log")
	parser.add_option("--os",action="store",type="string",dest="outStatistics",help="output file for DNA contamination statistics log")
	return parser

def mergeCutadaptLogs(cutadaptLog,outCutadapt):
    sample=os.path.split(cutadaptLog)[1].strip().split("_")[0]
    fout=open(outCutadapt,'w')
    with open(cutadaptLog,'r') as f:
        for line in f:
            if "Total reads processed" in line:
                total_reads=line.strip().split(":")[1].strip()
                fout.write(sample+"\t"+total_reads+"\t")
            elif "Reads with adapters" in line:
                reads_with_adapter=line.strip().split(":")[1].strip()
                fout.write(reads_with_adapter+"\t")
            elif "Reads that were too short" in line:
                reads_too_short=line.strip().split(":")[1].strip()
                fout.write(reads_too_short+"\t")
            elif "Reads written (passing filters)" in line:
                reads_left=line.strip().split(":")[1].strip()
                fout.write(reads_left+"\n")
            else:
                continue
    fout.close()

def mergeFilteringLogs(filteringLog,outFiltering):
    sample=os.path.split(filteringLog)[1].strip().split("_")[0]
    fout=open(outFiltering,'w')
    with open(filteringLog,'r') as f:
        for line in f:
            if "Input" in line:
                input_reads=line.strip().split(":")[1].strip().split(" ")[0]
                fout.write(sample+"\t"+input_reads+"\t")
            elif "Output" in line:
                Output_reads=line.strip().split(":")[1].strip().split(" ")[0]
                fout.write(Output_reads+"\t")
            elif "discarded" in line:
                discarded_reads=re.findall("\d+.*\(.*\)",line)[0]
                fout.write(discarded_reads+"\n")
            else:
                continue
    fout.close()

def mergerRNAContamLogs(rRNAcontamLog,outrRNAcontam):
    sample=os.path.split(rRNAcontamLog)[1].strip().split(".")[0]
    fout=open(outrRNAcontam,'w')
    with open(rRNAcontamLog,'r') as f:
        for line in f:
            if "reads processed" in line:
                processed_reads=line.strip().split(":")[1].strip()
                fout.write(sample+"\t"+processed_reads+"\t")
            elif "reads with at least one reported alignment" in line:
                reported_reads=line.strip().split(":")[1].strip()
                fout.write(reported_reads+"\t")
            elif "reads that failed to align" in line:
                failed_reads=line.strip().split(":")[1].strip()
                fout.write(failed_reads+"\n")
            else:
                continue
    fout.close()

def mergerMappingLogs(starMapping,outMapping):
    sample=os.path.split(starMapping)[1].strip().split(".")[0]
    fout=open(outMapping,'w')
    with open(starMapping,'r') as f:
        for line in f:
            if "Number of input reads" in line:
                input_reads=line.strip().split("|")[1].strip()
                fout.write(sample+"\t"+input_reads+"\t")
            elif "Uniquely mapped reads number" in line:
                unique_reads=line.strip().split("|")[1].strip()
            elif "Uniquely mapped reads %" in line:
                unique_ratio=line.strip().split("|")[1].strip()
                fout.write(unique_reads+" ("+unique_ratio+")"+"\t")
            elif "Number of reads mapped to too many loci" in line:
                multi_reads=line.strip().split("|")[1].strip()
            elif "% of reads mapped to too many loci" in line:
                multi_ratio=line.strip().split("|")[1].strip()
                fout.write(multi_reads+" ("+multi_ratio+")"+"\n")
            else:
                continue

    fout.close()

def mergerStatisticsLogs(statistics,outStatistics):
    sample=os.path.split(statistics)[1].strip().split("_")[0]
    fout=open(outStatistics,'w')
    with open(statistics,'r') as f:
        for line in f:
            if "unique mapped reads of exon" in line:
                exon_reads=line.strip().split(":")[1].strip()
                fout.write(sample+"\t"+exon_reads+"\t")
            elif "unique mapped reads of intergenic region" in line:
                DNA_reads=line.strip().split(":")[1].strip()
                fout.write(DNA_reads+"\t")
            elif "unique mapped reads of intron" in line:
                intron_reads=line.strip().split(":")[1].strip()
                fout.write(intron_reads+"\t")
            elif "unique mapped ambiguous reads of RNA" in line:
                RNA_ambiguous=line.strip().split(":")[1].strip()
                fout.write(RNA_ambiguous+"\n")
            else:
                continue
    fout.close()

def main():
    parser=create_parser_for_reads_summary()
    (options,args)=parser.parse_args()
    if options.cutadaptLog and options.outCutadapt:
        print("Start cutadapt...")
        mergeCutadaptLogs(options.cutadaptLog,options.outCutadapt)
    if options.filteringLog and options.outFiltering:
        print("Start filtering...")
        mergeFilteringLogs(options.filteringLog,options.outFiltering)
    if options.rRNAcontamLog and options.outrRNAcontam:
        print("Start rRNA contam...")
        mergerRNAContamLogs(options.rRNAcontamLog,options.outrRNAcontam)
    if options.starMapping and options.outMapping:
        print("Start star mapping...")
        mergerMappingLogs(options.starMapping,options.outMapping)
    if options.statistics and options.outStatistics:
        print("Start statistics...")
        mergerStatisticsLogs(options.statistics,options.outStatistics)
    print("Finish!")

if __name__=="__main__":
    main()