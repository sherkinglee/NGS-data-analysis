# ***Ribosome profiling（Ribo-seq）***

## **基本框架**

因为我从2018年4月开始定导在[清华大学杨雪瑞老师课题组](http://labyang.net/)之后接触的大部分项目都是基于Ribosome profiling数据分析来研究蛋白质的翻译调控的, 所以我的博士生阶段涉及到的主要研究数据类型就是ribosome profiling data，对该数据的分析流程也比较熟悉。这里我会将ribo-seq数据分析的大部分流程都介绍一遍，至于一些个性化的需求，则需要自己慢慢再探索。关于ribosome profling data analysis我会从以下几个方面进行介绍：

1）Ribosome profiling的简单背景介绍

2）Ribosome profiling数据类型和数据特点

3）Ribosome profiling数据上游处理流程

4）Ribosome profiling数据质控过程

5）Ribosome profiling数据下游处理过程


## **Ribosome profiling简单介绍**

[Ribosome profiling](https://en.wikipedia.org/wiki/Ribosome_profiling)也叫[Ribo-seq](https://en.wikipedia.org/wiki/Ribosome_profiling),是2009年[Ingolia等人](https://science.sciencemag.org/content/324/5924/218)结合深度测序技术，用来对活跃翻译的核糖体在mRNA分子上进行精确定位分析的工具。Ribosome profiling可以用于以下研究：

1）鉴定活跃翻译的开放式阅读框，即鉴定ORF；

2）鉴定翻译起始位点，寻找新的可能的翻译起始信号；

3）对蛋白质进行定量，一个polysome中可以结合多个核糖体，一个核糖体可以产生一条多肽链，用mRNA分子上结合的核糖体的数目表示某一时刻下该mRNA分子的翻译水平；可以用mRNA-seq做标准化，计算翻译效率（TE）表示单个mRNA分子上结合的核糖体数目；

4）监控翻译延长过程中出现的异常翻译调控事件，比如ribosome stalling，表示翻译延长速率收到抑制；

5）研究不同蛋白质之间的co-translation的情况。

## **Ribosome profiling数据类型和数据特点**

[Ribosome profiling](https://en.wikipedia.org/wiki/Ribosome_profiling)数据是基于Ribosome protected fragments测序达到的数据，所以原始数据类型也是fastq格式的。但是相比于mRNA-seq数据类型，ribo-seq数据具有以下几个特点

1）CDS富集： 因为是ribosome protected fragments,所以大部分都是CDS区域的reads

2）3-nt周期性：因为mRNA一般只会沿着一个固定的frame进行翻译，除非发生frame shift等特殊情况，所以在某个特定frame上的reads比其他两个frame上的要明显多很多，加上ribosome是没3个nt移动一个位置的，所以会呈现出明显的3-nt周期性；

3）start codon和stop codon出现峰值，主要是因为翻译起始速率和终止速率比较慢，当然也有CHX药物的作用；

4）长度分布： RNase I--28-30 nt, MNase---35nt

![image_1eermt8qu15bhqj8161q1a4vnl1m.png-197.5kB][1]

## **Ribosome profiling数据上游处理流程**
![image_1eern1mnv1ups1l8cnlk8vl12qk13.png-83.3kB][2]

### **00.rawdata, 数据下载**
```
## 00.rawdata, 数据下载
# GSE89704 as an example
for i in SRR50081{34..37};do fastq-dump $i;done
```
### **01.beforeQC,原始数据质控**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/01.beforeQC
fastaFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/00.rawdata

for i in SRR50081{34..37};
do
    fastqc $fastaFile/$i.fastq -o $workdir
done

```

### **02.cutadapt, 去3 end adapter**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/02.cutadapt
fastaFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/00.rawdata
adapter=CTGTAGGCACCATCAAT
#T_adapter=AGATCGGAAGAGCACACGTCT
#F_adapter=AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCG
for i in SRR50081{34..37}; ## RPF
do
        cutadapt -m 25 -M 45 --match-read-wildcards -a $adapter -o $workdir/$i.trimmed.fastq $fastaFile/$i.fastq > $workdir/${i}_trimmed.log

```

### **03.filter, remove reads with low quality**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/03.filter
fastaFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/02.cutadapt

for i in SRR50081{34..37};
do
    fastq_quality_filter -Q33 -v -q 25 -p 75 -i $fastaFile/$i.trimmed.fastq -o $workdir/$i.trimmed.Qfilter.fastq > $workdir/$i.Qfilter.log
done.
```

### **04.afterQC, filter后的质控**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200204/04.afterQC
fastaFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200204/03.filter

for i in SRR50081{34..37};
do
    fastqc $fastaFile/$i.trimmed.Qfilter.fastq -o $workdir
done

```

### **05.contam, remove rRNA contamination**

```
bowtie_noncoding_index=/Share2/home/lifj/Reference/human/hg38/rRNA_bowtieIndex/human_rRNA
fastqFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/03.filter
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/05.contam

for i in SRR50081{34..37};
do
    bsub -e $i.err -o $i.out -n 2 -q TEST-A "bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un=$workdir/noncontam_$i.fastq $bowtie_noncoding_index -q $fastqFile/$i.trimmed.Qfilter.fastq $workdir/$i.alin"
done

```

### **06.afterQC, final QC**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/06.finalQC
fastaFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/05.contam


for i in SRR50081{34..37};
do
    fastqc $fastaFile/noncontam_$i.fastq -o $workdir
done

```

### **07.STAR**

```
STAR_genome_index=/Share2/home/lifj/Reference/human/hg38/STAR_Human_Ensembl_GRCh38_Ensembl
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/07.STAR
fastqFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/05.contam


##  mapping
for i in SRR50081{34..37};
do
    mkdir -p $workdir/${i}_STAR
    STAR --runThreadN 8 --outFilterType Normal --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM --alignEndsType EndToEnd --out
done

## sort only transcriptome.out.bam
for i in SRR50081{34..37};
do
    samtools sort -T $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted -o $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.
done

## set index
for i in SRR50081{34..37};
do
    samtools index  $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam
    samtools index $workdir/${i}_STAR/$i.Aligned.sortedByCoord.out.bam
done

```

### **08.periodicity**

```
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/08.periodicity
BamDir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/07.STAR
RiboCode_annot=/Share2/home/lifj/Reference/human/hg38/RiboCode_annot_Human/RiboCode_annot

for i in SRR50081{34..37};
do
    metaplots -a $RiboCode_annot -r $BamDir/${i}_STAR/$i.Aligned.toTranscriptome.out.bam -o $workdir/$i -m 25 -M 45
done

```

### **09.readLengthDistribution**

```
pip install RiboMiner
workdir=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/09.readLengthDistribution
fastqFile=/Share2/home/lifj/Projects/04.XiaMen_DuanHaoran/RiboPipeline/20200709/03.filter

for i in SRR50081{34..37};
do
        LengthDistribution -i $fastqFile/$i.trimmed.Qfilter.fastq -o $workdir/$i
done
```

### **10. read count**

```
# pip install RiboMiner

ModifyHTseq -i bamFile -g gtfFile -o countsFile -t exon -m union -q 10 --minLen 25 --maxLen 35 --exclude-first 45 --exclude-last 15 --id-type gene_id

```

## **Ribosome profiling数据质控过程**

之前提到过Ribosome profiling data 有四个基本的特征，因此一个质量过关的Ribosome profiling数据也需要满足四个基本特征：

1） CDS富集
2） 3-nt周期性
3） start codon和stop codon有peak （非必须，但一般会有）
4）长度分布正常，RNase I消化得到的片段28-30nt，MNase消化的片段35nt

![image_1eesb4r7q5fas551kar1tgmgbq1g.png-62.6kB][3]


## **Ribosome profiling数据下游处理过程**

Ribosome profiling数据拿到mapping到基因组和转录组的bam文件之后，后续的分析有很多和传统RNA-seq数据很类似，包括利用ribosome profiling fragments (RPF)做蛋白质定量分析，用TE做翻译效率分析，即差异表达分析；除了这个之外，还可以用来鉴定ORF的数目和类别，尝试发现新的ORF;鉴定翻译起始位点，看看是否存在non-canonical翻译起始位点；还可以进行metagene 分析，看看是否会出现类似于Ribosome stalling这样的异常翻译调控事件；根据selective ribosome profiling数据，可以用来研究co-translation的翻译调控以及不同蛋白质之间的相互作用机制研究。

### **[差异表达分析](https://github.com/sherkinglee/NGS-data-analysis/blob/master/02.differential-expressioin-analysis/DEAnalysis.md)**

可以和RNA-seq数据一样，用DESeq2等差异表达软件来做差异表达分析，然后看看在某些条件下具体有哪些蛋白质的翻译发生了变化,详见[差异表达分析](https://github.com/sherkinglee/NGS-data-analysis/blob/master/02.differential-expressioin-analysis/DEAnalysis.md)专栏。但对于Ribo-seq数据来说，通常用翻译效率来做差异分析比较合适：详见[xtail](https://github.com/xryanglab/xtail/tree/master/vignettes).

```
## xtail 下载
    install.packages("devtools")
    library("devtools")
    install_github("xryanglab/xtail")
    
## xtail使用
library(xtail)

## 读入数据
mRNA <- read.table("mRNA.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
RPF <- read.table("RPF.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
common_gene <- intersect(rownames(mRNA),rownames(RPF))
mRNA <- mRNA[common_gene,]
RPF <- RPF[common_gene,]

## 构建condition 向量
condition <- rep(c("Normal","Tumor"),each=2)

## run xtail
TE_results <- xtail(mRNA,RPF,condition,bins=10000,minMeanCount = 10)
TE_DE_all <- resultsTable(TE_results)

```
拿到差异表达基因的list之后，可以做后续的功能富集分析以及其他类型的可视化分析，详见[火山图绘制](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Volcano_plots.md)，[热图绘制](https://github.com/sherkinglee/NGS-data-analysis/tree/master/05.plots)，[功能富集分析](https://github.com/sherkinglee/NGS-data-analysis/blob/master/02.differential-expressioin-analysis/visualization.md)等等。

### **ORF预测和翻译起始位点鉴定**

用Ribo-seq数据来做翻译起始位点鉴定和ORF预测的软件有很多，我们这里主要介绍我们实验室开发的[RiboCode](https://github.com/xryanglab/RiboCode)来完成ORF预测这一功能。

```
## 01.prepare annotation file
prepare_transcripts -g <Homo_sapiens.GRCh38.88.gtf> -f <Homo_sapiens.GRCh38.dna.primary_assembly.fa> -o <RiboCode_annot>

## 02. 3-nt periodicity
metaplots -a <RiboCode_annot> -r <transcript.bam> -o <output_prefix>

# 得到的configure文件格式
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-Ctrl-1   /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_ctrl_1_80S_final_STAR/si_ctrl_1_80S_final.Aligned.toTranscriptome.out.sorted.bam   yes 29,30,31    12,12,13

## 03.RiboCode
RiboCode -a <RiboCode_annot> -c <config.txt> -l no -g -o <RiboCode_ORFs_result>
```

![image_1eesdefug16aj24q1u2q16rm16rl1t.png-104.9kB][4]


**值得注意的是，用来绘制扇形图的数据文件是RiboCode_ORFs_result_collapsed.txt生成的，但是有些人会注意到该文件对应的ORF的数目和扇形图对应的数目不一样，原因是因为扇形图的数字表示含有ORF的基因的数目，而不是ORF的数目，所以只要把文件中的gene提取出来计算数目就会和扇形图的保持一致。另外，RPFcount计算出来的counts数有很多是0值，这个很奇怪，算是RiboCode的一个bug。**

上述是常规的使用方法，还有一种使用时将control和treat组的bam文件merge在一起，然后用来鉴定ORF,之后再统计control和treat组中这些鉴定出来的ORF的差异。

```
## 01. 构建new configure file
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-Ctrl-1   /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_ctrl_1_80S_final_STAR/si_ctrl_1_80S_final.Aligned.toTranscriptome.out.sorted.bam   yes 29,30,31    12,12,13
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-Ctrl-2   /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_ctrl_2_80S_final_STAR/si_ctrl_2_80S_final.Aligned.toTranscriptome.out.sorted.bam   yes 29,30,31    12,13,13
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-Ctrl-3   /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_ctrl_3_80S_final_STAR/si_ctrl_3_80S_final.Aligned.toTranscriptome.out.sorted.bam   yes 29,30,31    12,13,13
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-3e-1 /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_3e_1_80S_final_STAR/si_3e_1_80S_final.Aligned.toTranscriptome.out.sorted.bam   yes 29,30,31    12,13,13
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-3e-2 /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_3e_2_80S_final_STAR/si_3e_2_80S_final.Aligned.toTranscriptome.out.sorted.bam    yes 29,30,31    12,13,13
#SampleName     AlignmentFile   stranded(yes/reverse)   Psite   offset
si-3e-3 /workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/STARuniq/si_3e_3_80S_final_STAR/si_3e_3_80S_final.Aligned.toTranscriptome.out.sorted.bam    yes 29,30,31    12,13,13

## 02. run RiboCode
workdir=/workdata/LabMember2/lifj/lifj/Project/01.XiaMen_Human_EIF3E/From_Xingxd/Final_12samples/RiboCode/20191031
Ref=/workdata2/Projects/lifj/data/Reference/human/RiboCode_annot_Human/RiboCode_annot

RiboCode -a $Ref -c $workdir/config.txt -l yes -g -o RiboCode_ORFs_result -s ATG -A AAA,AAT,AAC,AAG,ATA,ATC,ATT,ACA,ACG,ACT,ACC,AGA,AGT,AGC,AGG,TAA,TAC,TAG,TAT,TTT,TTG,TTC,TTA,TCA,TCG,TCC,TCT,TGA,TGA,TGT,TGC,CAA,CAT,CAG,CAC,CTA,CTT,CTG,CTC,CCC,CCT,CCG,CCA,CGA,CGG,CGC,CGT,GAA,GAT,GAC,GAG,GTA,GTT,GTC,GTG,GCC,GCG,GCT,GCA,GGA,GGT,GGC,GGG -m 1
```

### **uORF的metagene plot**

有时候通过RiboCode预测出来了一些潜在的uORF，然后现在想看看uORF的周期性怎么样，可以用[RiboCode](https://github.com/xryanglab/RiboCode)自带的工具画图就可以，但是有时候想看看uORF上的 metagene  plot以及CDS区域的ORF情况。CDS的metaplot可以用[RiboMiner](https://github.com/xryanglab/RiboMiner)来实现，但是uORF的则可以通过[uORF_metagene.py](https://bitbucket.org/Sherking/ribopipe/src/master/RiboPipe/uORF_metagene.py)(私有库访问不了)来实现。还有一种需求是，计算5UTR区域的reads数据，然后看看和CDS区域counts的差异情况，可以用[UTR_counts.py](https://bitbucket.org/Sherking/ribopipe/src/master/RiboPipe/UTR_counts.py)(私有库)实现。


### **metagene analysis (图文不匹配，只是举个使用例子而已)**

Metagene Analysis可以用来发现那些潜在的翻译调控事件，这里可以用[RiboMiner](https://github.com/xryanglab/RiboMiner)来实现。RiboMiner具体的使用方法，详见[RiboMiner](https://github.com/xryanglab/RiboMiner)。
![image_1eesg0c60hsop601nhhg9jcoa34.png-247.7kB][5]

```
## 01. 构建configure.txt文件
bamFiles        readLengths     Offsets bamLegends
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_ctrl_1_80S_final_STAR/si_ctrl_1_80S_final.Aligned.toTranscriptome.out.sorted.bam     29,30,31        12,12,13        si_ctrl_1_80S
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_ctrl_2_80S_final_STAR/si_ctrl_2_80S_final.Aligned.toTranscriptome.out.sorted.bam     29,30,31        12,13,13        si_ctrl_2_80S
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_ctrl_3_80S_final_STAR/si_ctrl_3_80S_final.Aligned.toTranscriptome.out.sorted.bam     29,30,31        12,13,13        si_ctrl_3_80S
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_3e_1_80S_final_STAR/si_3e_1_80S_final.Aligned.toTranscriptome.out.sorted.bam 29,30,31        12,13,13        si_3e_1_80S
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_3e_2_80S_final_STAR/si_3e_2_80S_final.Aligned.toTranscriptome.out.sorted.bam 29,30,31        12,13,13        si_3e_2_80S
/Share2/home/lifj/Projects/03.RiboAnalyzer/XiaMen_test/STARuniq/si_3e_3_80S_final_STAR/si_3e_3_80S_final.Aligned.toTranscriptome.out.sorted.bam 29,30,31        12,13,13        si_3e_3_80S

## run metagene plot
MetageneAnalysis -f configure.txt -c ~/Reference/human/hg38/longest.transcripts.info.txt -o RocA03_up_CDS_unnormed -U codon -M RPKM -u 0 -d 500 -l 1 -n 1 -m 1 -e 5 --norm no -y 500 --CI 0.95 --type CDS -S 374_up_genes_info.txt --id-type transcript_id
```

![image_1eesg25do1c5ckai13em1km0uvi3h.png-475.7kB][6]

### **统计每个codon或者AA上的density，详见[RiboMiner](https://github.com/xryanglab/RiboMiner)**
### **计算不同转录本序列的序列特征，包括疏水性、带电性等等，详见[RiboMiner](https://github.com/xryanglab/RiboMiner)**
### **Co-translation研究，详见[RiboMiner](https://github.com/xryanglab/RiboMiner)**


  [1]: http://static.zybuluo.com/sherking/wjwkam9bzcdvy7rezuy1ggg8/image_1eermt8qu15bhqj8161q1a4vnl1m.png
  [2]: http://static.zybuluo.com/sherking/ngwdc0jemfz6atrb5yh8q5fz/image_1eern1mnv1ups1l8cnlk8vl12qk13.png
  [3]: http://static.zybuluo.com/sherking/n8vweyrht2us39n50zjbnvrh/image_1eesb4r7q5fas551kar1tgmgbq1g.png
  [4]: http://static.zybuluo.com/sherking/g7z0hhj2d1dfjtl3s9sxj6xq/image_1eesdefug16aj24q1u2q16rm16rl1t.png
  [5]: http://static.zybuluo.com/sherking/b69h9owi3tft3qg4hhnws5cb/image_1eesg0c60hsop601nhhg9jcoa34.png
  [6]: http://static.zybuluo.com/sherking/fzwksg4cbwy3rpg79j0j9ita/image_1eesg25do1c5ckai13em1km0uvi3h.png