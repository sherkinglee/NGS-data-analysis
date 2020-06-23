## **NGS数据前期处理过程**




## **Part 1: 数据下载**




目前常见的核酸数据库主要有三个，[NCBI](https://www.ncbi.nlm.nih.gov/)的[GeneBank](https://www.ncbi.nlm.nih.gov/)、日本的[DDBJ](https://www.ddbj.nig.ac.jp/index-e.html

)、还有[EMBL](https://www.embl.org/)的[EBI](https://www.ebi.ac.uk/)。最近，中国自己的核酸数据库[GSA](https://bigd.big.ac.cn/gsa/)。但是我只在NCBI上下载过原始数据，所以这里会用NCBI为例介绍数据的下载过程。




## **Download raw sequence data from  NCBI**




第一种方法：fastq-dump SRRID

首先下载sratoolkits

```
$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
$ tar zxvf sratoolkit.2.9.6-ubuntu64.tar.gz
$ cd sratoolkit.2.9.6-ubuntu64

#加入环境路径
$ echo 'export export PATH=$PATH:YOUR_PATH/sratoolkit.2.9.6-ubuntu64/bin' >> ~/.bash_profile
$ source ~/.bash_profile

```

然后进行数据的下载

```

## 这一步会直接将sra文件转换成fastq文件

fastq-dump SRR5515132

## 通过fastq-dump直接下载sra文件时会在~/ncbi/下产生中间sra文件，这需要手动将其删，如

workdir=/Share2/home/lifj/Projects/01.RiboPipe/human/GSE98623/00.rawdata

for i in SRR55151{32..49} ;do

    fastq-dump $i

    rm ~/ncbi/public/sra/$i.sra

don

## 对于PE数据来说，需要加一个参数--split-3 或者 --split-files

fastq-dump --split-3 SRR5515132

```

至于--split-files 和 --split-3有什么差别，网上有很多资料可以查到，比如说[here](https://www.jianshu.com/p/a8d70b66794c).

第二种方法：prefetch SRRID

```
## prefetch 下载下来的是sra文件，占用的空间更少，可以根据需要再转换成fastq，而且对于断点续传比较合适。

prefetch SRR5515132

## 也可以进行批量数据下载
## 将Sra ID号放进一个文件中，比如sra_list.txt
prefetch --option-file sra_list.txt

```

当然还存在其他的方法，比如说找到sra文件所在的目录，通过其URL利用wget进行下载就可以，网上也有很有这种介绍。


## **Download data from TCGA**

[TCGA](https://portal.gdc.cancer.gov/)项目的数据大家用的都比较多，很多时候我们也会从TCGA上下载一些数据，具体的下载方式也有很多，我们这里一般也是利用[GDC](https://portal.gdc.cancer.gov/)官方提供的方法进行下载。具体的下载方法见[GDC官网](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/).

第一步：下载gdc下载工具[gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)

![image_1e6lq4crf1reg12hc1a2k1e481rpcp.png-96.6kB][1]

第二步：根据[gdc-client的用法](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/)下载你感兴趣的文件即可。
```
## 根据UUID下载
gdc-client download 22a29915-6712-4f7a-8dba-985ae9a1f005

## 根据自己在DGC上挑选的文件manifest文件进行下载
dc-client download -m  /Users/JohnDoe/Downloads/gdc_manifest.txt

```

## **Part 2: 数据质控**

### **质控**
数据的质控过程，主要是看测序数据的质量，是否有adapter污染，测序深度怎样，reads长度多少，碱基比例如何，是否有rRNA等其它序列的污染等等。常用的测序数据质控软件主要有[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)、[multiqc](https://multiqc.info/)等等。

**FastQC使用**

```
## 第一步下载FastQC
略
## 第二步进行质控
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/01.beforeQC
fastaFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/00.rawdata

for i in SRR32088{52..57} SRR32088{64..69};
do
    fastqc $fastaFile/$i.fastq -o $workdir
done
```
FastQC的质控报告结果详见[FastQC官网](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

**[MultiQC](https://multiqc.info)使用**

```
## 第一步下载
略

## 第二步使用：multiqc是在做完所有前期数据分析之后才使用，只是一个汇总用的质控工具
multiqc input_dir -o output_dir

```
multiqc质控报告示例见[MultiQC官方网站](https://multiqc.info/examples/rna-seq/multiqc_report.html#)


## **Part 3: 数据清洗**

### **去除接头序列： cutadapt**

[cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)是比较常见的去除adapter的工具，可以用于双端也可以用于单端reads，但是一般用于单端会比较容易处理一些，PE双端数据我们可以用[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)或者[trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

第一步： cutadapt下载

```
## python pip下载
python3 -m pip install --user --upgrade cutadapt

## conda下载
conda create -n cutadaptenv cutadapt
```

第二步： cutadapt使用

```
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/02.cutadapt
fastaFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/00.rawdata
#adapter=CTGTAGGCACCATCAAT
adapter=AGATCGGAAGAGCACACGTCT

## try to remove the first 8 random nucleotides
for i in SRR32088{64..69}; ## RPF
do
        cutadapt -m 25 -M 35 --match-read-wildcards -a $adapter -o $workdir/$i.trimmed.fastq $fastaFile/$i.fastq > $workdir/${i}_trimmed.log

```

### **去除接头序列： Trimmomatic**

第一步： 下载[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```
下载到集群之后课将Trimmomatic的应用程序添加到PATH中，也可以使用绝对路径
```

第二步： [Trimmomatic使用](http://www.usadellab.org/cms/index.php?page=trimmomatic)

```
## pair end
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

## single end
java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### **去除接头序列： trim_galore**

第一步：下载[trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

```
略
```

第二步： 基本使用方法

```
workdir=/Share2/home/lifj/Projects/05.XieWei/20191112/02.cutadapt
rawdata=/Share2/home/lifj/Projects/05.XieWei/20191112/00.rawdata

for i in 1cell_Total_RNA_seq_rep1 2cell_early_72_total_RNA_seq_seq2 GV_100_total_RNA_seq late_2cell_Con_Total_RNA_seq_rep1 MII_100_total_RNA_seq_rep1;do
    mkdir -p $workdir/$i
    trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired $rawdata/${i}_r1.fq $rawdata/${i}_r2.fq --gzip -o $workdir/$i

done
```

### **去除低质量的reads： fastq_quality_filter**
**fastq_quality_filter**是[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)中的其中一个工具，用于过滤掉一些质量比较低的reads。

第一步： 下载[FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)
第二步：使用fastq_quality_filter

```
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/03.filter
fastaFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/02.cutadapt

for i in SRR32088{64..69} SRR32088{52..57};
do
    fastq_quality_filter -Q33 -v -q 25 -p 75 -i $fastaFile/$i.trimmed.fastq -o $workdir/$i.trimmed.Qfilter.fastq > $workdir/$i.Qfilter.log
done

```

### **去除rRNA等非编码RNA的污染**

第一步： 下载rRNA或者非编码RNA的核酸序列，可以从[ensemble](http://www.ensembl.org/info/data/ftp/index.html)中下载非编码RNA的fasta序列，也可以去[NCBI]()中搜索rRNA然后获取rRNA的序列。
第二步： 构建rRNA的bowtie或者bowtie2的index：

```
bowtie-build -f rRNA.fa human_rRNA
```

第三步： 将reads往rRNA上mapping

```
bowtie_noncoding_index=/workdata/LabMember2/lifj/lifj/data/Reference/human/rRNA_bowtieIndex/human_rRNA
fastqFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/03.filter
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE78959/05.contam

for i in SRR32088{64..69} SRR32088{52..57};
do
    bsub -e $i.err -o $i.out -n 2 "bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un=$workdir/noncontam_$i.fastq $bowtie_noncoding_index -q $fastqFile/$i.trimmed.Qfilter.fastq $workdir/$i.alin"
done

```

最后剩下来的reads可用于进一步的genome or transcriptome mapping的分析工作了。

  [1]: http://static.zybuluo.com/sherking/v0wy2e9v6rekwwvvkokpml3g/image_1e6lq4crf1reg12hc1a2k1e481rpcp.png