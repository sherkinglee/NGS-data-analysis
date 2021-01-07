# Snakemake实验流程搭建学习

标签：snakemake bioinformatics pipeline python3 20201220

[TOC]


## 概述

Snakemake是一个基于python3的流程搭建软件，主要是能够帮助我们搭建一些可以实现重复使用的分析pipeline。它的官方网址见[Snakemake](https://snakemake.readthedocs.io/en/stable/). 这是我第一次学习Snakemake，然后入门是基于[孟浩巍博士的视频课程](https://ke.qq.com/course/3031316?taid=10368983063478548)。大家感兴趣可以进去看一下。Snakemake有一些优点，可以方便我们进行常用生物信息流程的搭建，包括：

+ 基于python的流程搭建软件，语法友好；
+ 支持并行计算，不用手动写循环；
+ 支持断点运行，即使出错了，重新运行一下会自动跳过之前运行成功的部分，只运行断掉的部分；
+ 支持流程控制；
+ 支持CPU核心；
+ 支持运行时间控制；
+ more...

## [Snakemake安装](https://snakemake.readthedocs.io/en/stable/)

snakemake安装方式有很多，这里不再进行详细介绍，感兴趣的同学可以看官网[Snakemake安装](https://snakemake.readthedocs.io/en/stable/)。这里我使用conda进行安装：

```
## 先构建snakemake的虚拟环境
## first step: 安装anaconda3
## 添加镜像源
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
$ conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
# bioconda
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/

## second step: 创建虚拟环境，安装snakemake
$ conda create -n snakemake python=3.6 ## 创建snakemake的虚拟环境
$ conda activate snakemake ## 进入snakemake环境
$ conda install snakemake ## 安装snakemake
$ snakemake -h ## 测试是否安装成功
usage: snakemake [-h] [--dry-run] [--profile PROFILE]
                 [--cache [RULE [RULE ...]]] [--snakefile FILE] [--cores [N]]
                 [--local-cores N] [--resources [NAME=INT [NAME=INT ...]]]
                 [--set-threads RULE=THREADS [RULE=THREADS ...]]
                 [--set-scatter NAME=SCATTERITEMS [NAME=SCATTERITEMS ...]]
## 安装成功
$ conda deactivate ## 退出虚拟环境
```

## [Snakemake语法学习](https://snakemake.readthedocs.io/en/stable/)

snakemake基于不同的rule规则搭建一个分析流程，rule规则里面会定义好每个步骤的输入、输出、参数、执行命令等部分，将分析流程拆解为不同的组分。详细规则介绍见[Snakemake-tutorial](https://slides.com/johanneskoester/snakemake-tutorial).
<center>
<img src="http://static.zybuluo.com/sherking/nk4t8qki30yvems3cxew8a3p/image_1epvld40p1h3hknmphcj182tq9.png">
</center>

+ **mapping rule**

定义一个mapping的snakemake file，包含一个比对的rule， 名字为*bwa_map*.其中input包含两文件，一个*genome.fa* 以及一个测序文件*A.fastq*. 输出文件为bam文件，使用的命令为shell中定义的。
```
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

$ snakemake -np  ## -n 只判断逻辑对错，不执行； -p 输出执行命令

```

+ **sorting rule**

根据mapping的输出，定义一个sort rule, 利用samtools对上一步产生的文件进行sort。

```
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

+ **indel rule**

sort之后用samtools进行index的建立。

```
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```

完成比对之后呢，可以绘制拓扑图来观察上述过程的逻辑关系：

```
$ snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg
```

<center>
<img src="http://static.zybuluo.com/sherking/vadvag91v6b9cvtucpslxvo0/image_1epvn61lt1g221uto1t314qs1f6313.png">
</center>


## [搭建RiboMiner数据分析流程](https://github.com/xryanglab/RiboMiner)

先下载RiboMiner以及一些常用的二代测序数据软件：

+ [RiboMiner下载与安装](https://github.com/xryanglab/RiboMiner)

```
pip install RiboMiner
```

接下来，根据snakemake的语法规则，搭建一个RiboSeq数据分析平台。主要包括以下几个部分:

+ rule all (最后写）

rule all这个规则只有一个*input*，里面只定义最后需要保存的文件或者目录，用户可以根据需要选择保留那些文件。

```
ruleorder:statistic_contamination>summary>mergeLogs
rule all:
    input:
        expand("01.beforeQC/{sample}",sample=SAMPLES),
        expand("02.cutadapt/{sample}_trimmed.fastq",sample=SAMPLES),
        expand("03.filter/{sample}_trimmed_Qfilter.fastq",sample=SAMPLES),
        expand("04.afterQC/{sample}",sample=SAMPLES),
        expand("05.contam/noncontam_{sample}.fastq",sample=SAMPLES),
        expand("06.finalQC/{sample}",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Log.final.out",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai",sample=SAMPLES,type={"sortedByCoord.out","toTranscriptome.out.sorted"}),
        expand("07.STAR/{sample}_STAR/{sample}.bw",sample=SAMPLES),
        expand("08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",sample=SAMPLES),
        expand("08.periodicity/{sample}_pre_config.txt",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.pdf",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.txt",sample=SAMPLES),
        expand("10.read_counts/{sample}.counts",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_DNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_Intron.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_RNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_reads_distribution.txt",sample=SAMPLES),
        expand("{sample}_Cutadapt.txt",sample=SAMPLES),
        expand("{sample}_Filtering.txt",sample=SAMPLES),
        expand("{sample}_RRNA_contam.txt",sample=SAMPLES),
        expand("{sample}_Star_mapping.txt",sample=SAMPLES),
        expand("{sample}_Statistics.txt",sample=SAMPLES),
        expand("12.summary/Ribo_seq_summary_statistics.txt")
```

+ samples, paramters and tools defination （根据流程需要随时写）

```
import os

## samples defination
SAMPLES=["SRR8445000",
         "SRR8445001",
	 "SRR8445002",
	 "SRR8445003"]


## parameters defination
BOWTIE_NONCODING_INDEX="/workdata/home/lifj/lifj/data/Reference/human/rRNA_bowtieIndex/human_rRNA"
STAR_GENOME_INDEX="/workdata/home/lifj/lifj/data/Reference/human/STAR_Human_Ensembl_GRCh38_Ensembl"
RIBOCODE_ANNOT="/workdata/home/lifj/lifj/data/Reference/human/RiboCode_annot_Human/RiboCode_annot"
GTFFILE="/workdata/home/lifj/lifj/data/Reference/human/Homo_sapiens.GRCh38.88.gtf"
ADAPTER="TGGAATTCTCGGGTGCCA"


## softwares defination
# make sure the all softwares needed are all installed and in your PATH

with os.popen("which fastqc") as path:
    FASTQC = path.read().strip()
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
with os.popen("which fastq_quality_filter") as path:
    FILTER = path.read().strip()
with os.popen("which bowtie") as path:
    BOWTIE = path.read().strip()
with os.popen("which STAR") as path:
    STAR = path.read().strip()
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
with os.popen("which metaplots") as path:
    METAPLOTS = path.read().strip()
with os.popen("which LengthDistribution") as path:
    LENGTHDISTRIBUTION = path.read().strip()
with os.popen("which ModifyHTseq") as path:
    MODIFYHTSEQ = path.read().strip()
with os.popen("which StatisticReadsOnDNAsContam") as path:
    STATISTICREADSONDNASCONTAM = path.read().strip()
with os.popen("which bamCoverage") as path:
    BAMCOVERAGE = path.read().strip()

```

+ beforeQC

```
rule beforeQC:
    input:
        "../00.rawdata/{sample}.fastq"
    output:
        directory("01.beforeQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """
```

+ cutadapt

```
rule Cutadapt:
    input:
        "../00.rawdata/{sample}.fastq"
    output:
        "02.cutadapt/{sample}_trimmed.fastq"
    log:
        "02.cutadapt/{sample}_trimmed.log"
    threads: 1
    shell:
        """
        {CUTADAPT} -m 20 -M 40 --match-read-wildcards -a {ADAPTER} -o {output} {input} > {log} 2>&1
        """
```

+ quality filter

```
rule quality_filter:
    input:
        "02.cutadapt/{sample}_trimmed.fastq"
    output:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    log:
        "03.filter/{sample}_trimmed_Qfilter.log"
    threads: 1
    shell:
        """
        {FILTER} -Q33 -v -q 25 -p 75 -i {input} -o {output} > {log} 2>&1
        """
```

+ afterQC

```
rule afterQC:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        directory("04.afterQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """
## 若出现Can't exec "java": Permission denied at /Share2/home/lifj/Software/FastQC/fastqc line 283
## 这是个很迷的报错，我即使修改了fastqc的权限，依然报错，不知道是哪出现问题。
```

+ remove_rRNA_contamination

```
rule remove_rRNA:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        "05.contam/{sample}.aln",
        "05.contam/noncontam_{sample}.fastq"
    log:
        "05.contam/{sample}.log"
    threads: 1
    shell:
        """
        {BOWTIE} -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un={output[1]} {BOWTIE_NONCODING_INDEX} -q {input} {output[0]} > {log} 2>&1
        """
```

+ finalQC

```
rule finalQC:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        directory("06.finalQC/{sample}")
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """
```

+ STAR mapping

```
rule STAR:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Log.final.out",
        "07.STAR/{sample}_STAR/{sample}.Log.out",
        "07.STAR/{sample}_STAR/{sample}.Log.progress.out",
        "07.STAR/{sample}_STAR/{sample}.ReadsPerGene.out.tab",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.SJ.out.tab"
    threads: 4
    shell:
        """
        {STAR} --runThreadN {threads} --outFilterType Normal --outWigType wiggle \
        --outWigStrand Stranded --outWigNorm RPM --alignEndsType EndToEnd \
        --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir {STAR_GENOME_INDEX} \
        --readFilesIn {input} \
        --outFileNamePrefix 07.STAR/{wildcards.sample}_STAR/{wildcards.sample}. \
        --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
        --outSAMattributes All
        """


## 若是出现unique mapping reads特别低的并且star log file中显示Number of reads unmapped: too short太多时，可以调节--outFilterMatchNmin 参数，如换做15 or 16那么这部分比例会降低，unqiue和multiple mapping的reads会升高。
```



+ samtools sort

```
rule samtools_sort:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    threads: 1
    shell:
        "{SAMTOOLS} sort -T {output}.temp -o {output} {input}"
```

+ samtools index

```
rule samtools_index:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai"
    threads: 1
    shell:
        "{SAMTOOLS} index {input}"
```

+ bam2bigwig

```
rule bam2bigwig:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.bw"
    log:
        "07.STAR/{sample}_STAR/{sample}.bam2bw.log"
    threads: 2
    shell:
        """
        {BAMCOVERAGE} -b {input} -o {output} --normalizeUsing CPM --binSize 1 > {log} 2>&1
        """
## 这一步比较耗时
```

+ 3-nt periodicity checking

```
rule periodicity:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",
        "08.periodicity/{sample}_pre_config.txt"
    threads: 1
    shell:
        "{METAPLOTS} -a {RIBOCODE_ANNOT} -r {input} -o 08.periodicity/{wildcards.sample} -m 20 -M 40"
```

+ length statistics

```
rule length_distribution:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    output:
        "09.LengthDistribution/{sample}_reads_length.pdf",
        "09.LengthDistribution/{sample}_reads_length.txt"
    log:
        "09.LengthDistribution/{sample}_reads_length.log"
    threads: 1
    shell:
        "{LENGTHDISTRIBUTION} -i {input} -o 09.LengthDistribution/{wildcards.sample} -f bam > {log} 2>&1"

```

+ htseq_counts

```
rule htseq_count:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "10.read_counts/{sample}.counts"
    log:
        "10.read_counts/{sample}.counts.log"
    threads: 1
    shell:
        """
        {MODIFYHTSEQ} -i {input} -g {GTFFILE} -o {output} -t CDS -m union -q 10 --minLen 25 --maxLen 35 --exclude-first 45 --exclude-last 15 --id-type gene_name > {log} 2>&1
        """

```

+ statistic contaminations

```
rule statistic_contamination:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "11.RPF_statistics/{sample}_DNA.pdf",
        "11.RPF_statistics/{sample}_Intron.pdf",
        "11.RPF_statistics/{sample}_RNA.pdf",
        "11.RPF_statistics/{sample}_reads_distribution.txt"
    log:
        "11.RPF_statistics/{sample}.statistics.log"
    threads: 1
    shell:
        """
        {STATISTICREADSONDNASCONTAM} -i  {input}  -g {GTFFILE}  -o  11.RPF_statistics/{wildcards.sample} > {log} 2>&1
        """
```


+ summary log files of each step

```
rule summary:
    input:
        cutadapt="02.cutadapt/{sample}_trimmed.log",
        filtering="03.filter/{sample}_trimmed_Qfilter.log",
        removeRNAcontam="05.contam/{sample}.log",
        starMapping="07.STAR/{sample}_STAR/{sample}.Log.final.out",
        statistics="11.RPF_statistics/{sample}_reads_distribution.txt"
    output:
        cutadaptLog=temp("{sample}_Cutadapt.txt"),
        filteringLog=temp("{sample}_Filtering.txt"),
        rRNAContamLog=temp("{sample}_RRNA_contam.txt"),
        starMappingLog=temp("{sample}_Star_mapping.txt"),
        statisticsLog=temp("{sample}_Statistics.txt")
        # finalLog="RNA_seq_summary_statistics.txt"

    shell:
        """
        python summary.py -f {input.filtering}  --of {output.filteringLog} \
        -c {input.cutadapt}  --oc {output.cutadaptLog} \
        -r {input.removeRNAcontam} --or {output.rRNAContamLog} \
        -s {input.statistics} --os {output.statisticsLog}\
        -m {input.starMapping} --om {output.starMappingLog}
        """


## 其中summary.py的内如如下：
$ vim summary.py
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

```

+ merge all logs files

```
rule mergeLogs:
    input:
        cutadaptLogs=[filename for filename in os.listdir(".") if "Cutadapt" in  filename],
        filteringLogs=[filename for filename in os.listdir(".") if "Filtering" in  filename],
        rRNAContamLogs=[filename for filename in os.listdir(".") if "RRNA_contam" in  filename],
        starMappingLogs=[filename for filename in os.listdir(".") if "Star_mapping" in  filename],
        statisticsLogs=[filename for filename in os.listdir(".") if "Statistics" in  filename]
    output:
        "12.summary/RNA_seq_summary_statistics.txt"
    run:
        shell("echo -e '# cutadapt\nsample\tTotal\tTrimmed(Percent)\tshortNum(Percentage)\tLeftNum(Percentage)' >> {output} ")
        shell("cat {input.cutadaptLogs}>> {output}")
        shell("echo -e '# filtering\nsample\tinputNum\tRemained\tdiscarded(Percent)' >> {output} ")
        shell("cat {input.filteringLogs}>> {output}")
        shell("echo -e '# remove RNA contamination\nsample\tProcessedNum\trRNA(Percent)\tnoContamRNA(Percent)' >> {output} ")
        shell("cat {input.rRNAContamLogs}>> {output}")
        shell("echo -e '# Star mapping\nsample\tinput\tUniquelyMapped(Percent)\tMutipulMapped(Percent)' >> {output} ")
        shell("cat {input.starMappingLogs}>> {output}")
        shell("echo -e '# DNA contamination\nsample\tExon\tDNA\tIntron\tambiguous_RNA' >> {output} ")
        shell("cat {input.statisticsLogs}>> {output}")

```


+ summary all rules of snakemake file

最后整个snakemake流程如下：

```
# -*- coding:UTF-8 -*-
'''
Author: Li Fajin
Date: 2020-12-21 19:59:09
LastEditors: Li Fajin
LastEditTime: 2021-01-01 16:56:30
Description: snakemake pipeline for ribo-seq data analyses, just for basic control steps such fastqc, cutadapt, qfiltering, mapping, samtools sort and index, 3nt periodicity checking, et al.
'''

import os

## samples defination
SAMPLES=["SRR8445000",
         "SRR8445001",
	 "SRR8445002",
	 "SRR8445003"]


## parameters defination
BOWTIE_NONCODING_INDEX="/workdata/home/lifj/lifj/data/Reference/human/rRNA_bowtieIndex/human_rRNA"
STAR_GENOME_INDEX="/workdata/home/lifj/lifj/data/Reference/human/STAR_Human_Ensembl_GRCh38_Ensembl"
RIBOCODE_ANNOT="/workdata/home/lifj/lifj/data/Reference/human/RiboCode_annot_Human/RiboCode_annot"
GTFFILE="/workdata/home/lifj/lifj/data/Reference/human/Homo_sapiens.GRCh38.88.gtf"
ADAPTER="TGGAATTCTCGGGTGCCA"


## softwares defination
# make sure the all softwares needed are all installed and in your PATH

with os.popen("which fastqc") as path:
    FASTQC = path.read().strip()
with os.popen("which cutadapt") as path:
    CUTADAPT = path.read().strip()
with os.popen("which fastq_quality_filter") as path:
    FILTER = path.read().strip()
with os.popen("which bowtie") as path:
    BOWTIE = path.read().strip()
with os.popen("which STAR") as path:
    STAR = path.read().strip()
with os.popen("which samtools") as path:
    SAMTOOLS = path.read().strip()
with os.popen("which metaplots") as path:
    METAPLOTS = path.read().strip()
with os.popen("which LengthDistribution") as path:
    LENGTHDISTRIBUTION = path.read().strip()
with os.popen("which ModifyHTseq") as path:
    MODIFYHTSEQ = path.read().strip()
with os.popen("which StatisticReadsOnDNAsContam") as path:
    STATISTICREADSONDNASCONTAM = path.read().strip()
with os.popen("which bamCoverage") as path:
    BAMCOVERAGE = path.read().strip()


## snakemakes pipeline
ruleorder:statistic_contamination>summary>mergeLogs ## 如果不加这一个，整个流程运行下来会报错，因为会优先执行mergeLogs
rule all:
    input:
        expand("01.beforeQC/{sample}",sample=SAMPLES),
        expand("02.cutadapt/{sample}_trimmed.fastq",sample=SAMPLES),
        expand("03.filter/{sample}_trimmed_Qfilter.fastq",sample=SAMPLES),
        expand("04.afterQC/{sample}",sample=SAMPLES),
        expand("05.contam/noncontam_{sample}.fastq",sample=SAMPLES),
        expand("06.finalQC/{sample}",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Log.final.out",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam",sample=SAMPLES),
        expand("07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai",sample=SAMPLES,type={"sortedByCoord.out","toTranscriptome.out.sorted"}),
        expand("07.STAR/{sample}_STAR/{sample}.bw",sample=SAMPLES),
        expand("08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",sample=SAMPLES),
        expand("08.periodicity/{sample}_pre_config.txt",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.pdf",sample=SAMPLES),
        expand("09.LengthDistribution/{sample}_reads_length.txt",sample=SAMPLES),
        expand("10.read_counts/{sample}.counts",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_DNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_Intron.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_RNA.pdf",sample=SAMPLES),
        expand("11.RPF_statistics/{sample}_reads_distribution.txt",sample=SAMPLES),
        expand("{sample}_Cutadapt.txt",sample=SAMPLES),
        expand("{sample}_Filtering.txt",sample=SAMPLES),
        expand("{sample}_RRNA_contam.txt",sample=SAMPLES),
        expand("{sample}_Star_mapping.txt",sample=SAMPLES),
        expand("{sample}_Statistics.txt",sample=SAMPLES),
        expand("12.summary/Ribo_seq_summary_statistics.txt")

rule beforeQC:
    input:
        "../00.rawdata/{sample}.fastq"
    output:
        directory("01.beforeQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule Cutadapt:
    input:
        "../00.rawdata/{sample}.fastq"
    output:
        "02.cutadapt/{sample}_trimmed.fastq"
    log:
        "02.cutadapt/{sample}_trimmed.log"
    threads: 1
    shell:
        """
        {CUTADAPT} -m 20 -M 40 --match-read-wildcards -a {ADAPTER} -o {output} {input} > {log} 2>&1
        """

rule quality_filter:
    input:
        "02.cutadapt/{sample}_trimmed.fastq"
    output:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    log:
        "03.filter/{sample}_trimmed_Qfilter.log"
    threads: 1
    shell:
        """
        {FILTER} -Q33 -v -q 25 -p 75 -i {input} -o {output} > {log} 2>&1
        """

rule afterQC:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        directory("04.afterQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule remove_rRNA:
    input:
        "03.filter/{sample}_trimmed_Qfilter.fastq"
    output:
        "05.contam/{sample}.aln",
        "05.contam/noncontam_{sample}.fastq"
    log:
        "05.contam/{sample}.log"
    threads: 1
    shell:
        """
        {BOWTIE} -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un={output[1]} {BOWTIE_NONCODING_INDEX} -q {input} {output[0]} > {log} 2>&1
        """


rule finalQC:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        directory("06.finalQC/{sample}")
    threads: 1
    shell:
        """
        mkdir -p {output}
        {FASTQC} {input} -o {output}
        """

rule STAR:
    input:
        "05.contam/noncontam_{sample}.fastq"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam",
        "07.STAR/{sample}_STAR/{sample}.Log.final.out",
        "07.STAR/{sample}_STAR/{sample}.Log.out",
        "07.STAR/{sample}_STAR/{sample}.Log.progress.out",
        "07.STAR/{sample}_STAR/{sample}.ReadsPerGene.out.tab",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.UniqueMultiple.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str1.out.wig",
        "07.STAR/{sample}_STAR/{sample}.Signal.Unique.str2.out.wig",
        "07.STAR/{sample}_STAR/{sample}.SJ.out.tab"
    threads: 4
    shell:
        """
        ## part of parameters from GSE125114
        {STAR} --runThreadN {threads} --outFilterType Normal --outWigType wiggle --outWigStrand Stranded \
            --outWigNorm RPM --alignEndsType EndToEnd --outFilterMismatchNmax 2 \
            --outFilterMultimapNmax 20 --outFilterMismatchNoverLmax 0.04 --outFilterMatchNmin 15 --genomeDir {STAR_GENOME_INDEX} --readFilesIn {input} \
            --outFileNamePrefix  07.STAR/{wildcards.sample}_STAR/{wildcards.sample}. \
            --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts \
            --outSAMattributes All
        """

rule samtools_sort:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    threads: 1
    shell:
        "{SAMTOOLS} sort -T {output}.temp -o {output} {input}"


rule samtools_index:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.Aligned.{type}.bam.bai"
    threads: 1
    shell:
        "{SAMTOOLS} index {input}"

rule bam2bigwig:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "07.STAR/{sample}_STAR/{sample}.bw"
    log:
        "07.STAR/{sample}_STAR/{sample}.bam2bw.log"
    threads: 2
    shell:
        """
        {BAMCOVERAGE} -b {input} -o {output} --normalizeUsing CPM --binSize 1 > {log} 2>&1
        """

rule periodicity:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.bam"
    output:
        "08.periodicity/{sample}{sample}.Aligned.toTranscriptome.out.pdf",
        "08.periodicity/{sample}_pre_config.txt"
    threads: 1
    shell:
        "{METAPLOTS} -a {RIBOCODE_ANNOT} -r {input} -o 08.periodicity/{wildcards.sample} -m 20 -M 40"

rule length_distribution:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.toTranscriptome.out.sorted.bam"
    output:
        "09.LengthDistribution/{sample}_reads_length.pdf",
        "09.LengthDistribution/{sample}_reads_length.txt"
    log:
        "09.LengthDistribution/{sample}_reads_length.log"
    threads: 1
    shell:
        "{LENGTHDISTRIBUTION} -i {input} -o 09.LengthDistribution/{wildcards.sample} -f bam > {log} 2>&1"


# specifically for Ribo-seq data, if you have RNA-seq data, using htseq-count please.
rule htseq_count:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "10.read_counts/{sample}.counts"
    log:
        "10.read_counts/{sample}.counts.log"
    threads: 1
    shell:
        """
        {MODIFYHTSEQ} -i {input} -g {GTFFILE} -o {output} -t CDS -m union -q 10 --minLen 25 --maxLen 35 --exclude-first 45 --exclude-last 15 --id-type gene_name > {log} 2>&1
        """

rule statistic_contamination:
    input:
        "07.STAR/{sample}_STAR/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        "11.RPF_statistics/{sample}_DNA.pdf",
        "11.RPF_statistics/{sample}_Intron.pdf",
        "11.RPF_statistics/{sample}_RNA.pdf",
        "11.RPF_statistics/{sample}_reads_distribution.txt"
    log:
        "11.RPF_statistics/{sample}.statistics.log"
    threads: 1
    shell:
        """
        {STATISTICREADSONDNASCONTAM} -i  {input}  -g {GTFFILE}  -o  11.RPF_statistics/{wildcards.sample} > {log} 2>&1
        """


rule summary:
    input:
        cutadapt="02.cutadapt/{sample}_trimmed.log",
        filtering="03.filter/{sample}_trimmed_Qfilter.log",
        removeRNAcontam="05.contam/{sample}.log",
        starMapping="07.STAR/{sample}_STAR/{sample}.Log.final.out",
        statistics="11.RPF_statistics/{sample}_reads_distribution.txt"
    output:
        cutadaptLog=temp("{sample}_Cutadapt.txt"),
        filteringLog=temp("{sample}_Filtering.txt"),
        rRNAContamLog=temp("{sample}_RRNA_contam.txt"),
        starMappingLog=temp("{sample}_Star_mapping.txt"),
        statisticsLog=temp("{sample}_Statistics.txt")
        # finalLog="RNA_seq_summary_statistics.txt"

    shell:
        """
        python summary.py -f {input.filtering}  --of {output.filteringLog} -c {input.cutadapt}  --oc {output.cutadaptLog} \
            -r {input.removeRNAcontam} --or {output.rRNAContamLog} -s {input.statistics} --os {output.statisticsLog}\
            -m {input.starMapping} --om {output.starMappingLog}
        """

rule mergeLogs:
    input:
        cutadaptLogs=[filename for filename in os.listdir(".") if "Cutadapt" in  filename],
        filteringLogs=[filename for filename in os.listdir(".") if "Filtering" in  filename],
        rRNAContamLogs=[filename for filename in os.listdir(".") if "RRNA_contam" in  filename],
        starMappingLogs=[filename for filename in os.listdir(".") if "Star_mapping" in  filename],
        statisticsLogs=[filename for filename in os.listdir(".") if "Statistics" in  filename]
    output:
        "12.summary/RNA_seq_summary_statistics.txt"
    run:
        shell("echo -e '# cutadapt\nsample\tTotal\tTrimmed(Percent)\tshortNum(Percentage)\tLeftNum(Percentage)' >> {output} ")
        shell("cat {input.cutadaptLogs}>> {output}")
        shell("echo -e '# filtering\nsample\tinputNum\tRemained\tdiscarded(Percent)' >> {output} ")
        shell("cat {input.filteringLogs}>> {output}")
        shell("echo -e '# remove RNA contamination\nsample\tProcessedNum\trRNA(Percent)\tnoContamRNA(Percent)' >> {output} ")
        shell("cat {input.rRNAContamLogs}>> {output}")
        shell("echo -e '# Star mapping\nsample\tinput\tUniquelyMapped(Percent)\tMutipulMapped(Percent)' >> {output} ")
        shell("cat {input.starMappingLogs}>> {output}")
        shell("echo -e '# DNA contamination\nsample\tExon\tDNA\tIntron\tambiguous_RNA' >> {output} ")
        shell("cat {input.statisticsLogs}>> {output}")



```

## 运行snakemake file

如何运行这个snakefile文件呢，可以在本地运行，也可以提交到服务器运行，因为涉及到mapping等耗时过程，所以需要提交到服务器进行运行。

+ 测试逻辑是否正确

```
## -n dry-run，测试逻辑，不运行； -p打印每步的命令;这步会输出任务数目
$snakemake -s riboseq-snakemake.py -np
......
Job counts:
        count   jobs
        2       Cutadapt
        2       STAR
        2       afterQC
        1       all
        2       beforeQC
        2       finalQC
        2       htseq_count
        2       length_distribution
        2       periodicity
        2       quality_filter
        2       remove_rRNA
        4       samtools_index
        2       samtools_sort
        2       bam2bigwig
        2       statistic_contamination
        31
This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

```

+ 本地运行 （不推荐）

```
$ snakemake -s riboseq-snakemake.py -j 8 ## 核数
$ snakemake -s riboseq-snakemake.py -j 8 --rerun-incomplete
```

+ 提交到计算节点运行（LSF提交系统， bsub提交）

```
## 先开一个screen
$ screen -S riboseq
$ snakemake -s riboseq-snakemake.py --cluster "bsub -q TEST-A -n 4 -e err -o out" --jobs 31

##或者编写提交脚本
$ vim riboseq_bsub.sh
#########################################################################
# File Name: riboseq_bsub.sh
# Author: lifj
# mail: lfj17@mails.tsinghua.edu.cn
# Created Time: Mon 21 Dec 2020 11:06:05 AM CST
#########################################################################
#!/bin/bash
#BSUB -J riboseq_bsub.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 4
#BSUB -q TEST-A

snakemake -s riboseq-snakemake.py --cluster "bsub -q TEST-A -o err -e out -n 4" --jobs 31

$ bsub < riboseq_bsub.sh
## 一般不会有问题，但是有些任务会因为不知名原因断掉，如fastqc那一步，这个可以在整个pipeline结束之后再通过以下运行一次，可以跑完

$ snakemake -s riboseq-snakemake.py --cluster "bsub -q TEST-A -o err -e out -n 4" --jobs 10 #根据需要改写

```

+ 提交到集群运行snakemake (slurm提交系统， sbatch提交）

```
## 先开一个screen
$ screen -S a
$ snakemake -s Ribo-seq-snakemake.py  --cluster "sbatch -p cmp -N 2 -n 8 -e Ribo.err -o Ribo.out " --jobs 61

```

<center>
<img src="http://static.zybuluo.com/sherking/valorut9817sqjbpv68m9mxb/image_1epvhfmve1g3inbm18d51so9lm49.png">
</center>

+ 流程拓扑图

```
$ snakemake -s riboseq-snakemake.py --dag | dot -Tpdf > riboseq.pdf
```

<center>
<img src="http://static.zybuluo.com/sherking/r7s45ud2o5f98dd6700l1c0a/image_1eq4rnnq8182e10eu1bbu8utj8u9.png">
</center>
