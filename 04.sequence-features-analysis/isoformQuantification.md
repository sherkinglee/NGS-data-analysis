## **转录本isoform定量**

常见的转录本定量工具有很多，包括所谓的alignment-based，如RSEM等，也包括alignment-free的如Salmon等等。在研究计算K562细胞中ZAK两个isoform表达量的时候可能需要用Salmon测试一下，所以这里先记录一下Salmon进行转录本定量的使用方法。

### **[Salmon进行转录本定量的原理](http://robpatro.com/redesign/Quantification.pdf)**

或者参见[Salmon原理](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/04_quasi_alignment_salmon.html)

![image_1eg2om6bim3s2ude3dlk33eg9.png-190.5kB][1]
![image_1eg2omu413j3tf0c53163oo7hm.png-207.7kB][2]


### **[Salmon安装](https://combine-lab.github.io/salmon/getting_started/)**

可以根据[Salmon官网](https://combine-lab.github.io/salmon/getting_started/)进行安装，我是直接用conda直接安装的

```
conda install -c bioconda salmon
```

尽管顺利安装上了，但是在使用的时候发现并没有安装成功，因为出现如下错误：

```
$ salmon -h
salmon: error while loading shared libraries: libboost_iostreams.so.1.60.0: cannot open shared object file: No such file or directory

```
通过以下命令发现，libboost_isostreams的版本太老了：

```
$ locate libboost_iostreams
/opt/UCSF/Chimera64-1.10.1/lib/libboost_iostreams.so
/opt/UCSF/Chimera64-1.10.1/lib/libboost_iostreams.so.1.52.0
```

考虑到自己没有权限重新安装新版本的package，所以直接从[github](https://github.com/COMBINE-lab/salmon/releases/tag/v1.3.0)中下载tar包进行安装：

```
$ wget -c https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz

$ tar -xvzf salmon-1.3.0_linux_x86_64.tar.gz
$ cd salmon-latest_linux_x86_64/bin
$ ./salmon -h
salmon v1.3.0

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     index      : create a salmon index
     quant      : quantify a sample
     alevin     : single cell analysis
     swim       : perform super-secret operation
     quantmerge : merge multiple quantifications into a single file
```

### **[Salmon使用](https://combine-lab.github.io/salmon/getting_started/)**

#### **1.下载转录组序列**

原本我以为需要用StringTie或者trinity等软件进行组装，但是根据salmon官网的介绍，发现可以从[Ensemble geome browser](https://asia.ensembl.org/info/data/ftp/index.html)直接下载每个物种的cDNA序列即可，以人hg38为例。

```
$ wget -c  ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

#### **2.为转录组建index**
![image_1eg2qpqop2am1agc18h611c9if813.png-34.9kB][3]

```
salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i hg38 -k 31 ## PE >=75 nt 最优，如果reads短，可以将-k缩小

salmon index -t Homo_sapiens.GRCh38.cdna.all.fa.gz -i hg38 -k 21 ## 50nt
```

#### **3.下载RNA-seq数据进行测试**

这里下载[ENCODE](https://www.encodeproject.org/)中K562细胞的PE RNA-seq数据进行测试。

```
$ cat UTL.txt
https://www.encodeproject.org/files/ENCFF549JWM/@@download/ENCFF549JWM.fastq.gz
https://www.encodeproject.org/files/ENCFF923HGD/@@download/ENCFF923HGD.fastq.gz
https://www.encodeproject.org/files/ENCFF154VOX/@@download/ENCFF154VOX.fastq.gz
https://www.encodeproject.org/files/ENCFF119MRB/@@download/ENCFF119MRB.fastq.gz

$ cat download.sh
cat URL.txt|while read line
do
wget  -c $line
done

$ bash download.sh

```

#### **4.转录本计数**

若是单端数据：

```
$ salmon quant -i hg38 -l A -r sample.fastq -o outputdir --useVBOpt --seqBias --validateMappings

```

我这里用的是PE-50的双端数据：

```
#!/bin/bash
#BSUB -J run.sh
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 4
#BSUB -q TEST-A

workdir=/Share2/home/lifj/Projects/11.ENCODE/fastq/K562/02.Salmon
fastq=/Share2/home/lifj/Projects/11.ENCODE/fastq/K562/00.rawdata
index=/Share2/home/lifj/Reference/human/hg38/Salmon_index/hg38-k21

salmon quant -i $index -l A -1 $fastq/ENCFF549JWM.fastq -2 $fastq/ENCFF923HGD.fastq --useVBOpt --seqBias --validateMappings -o NWK_transcripts_quant

salmon quant -i $index -l A -1 $fastq/ENCFF154VOX.fastq -2 $fastq/ENCFF119MRB.fastq --useVBOpt --seqBias --validateMappings -o NSR_transcripts_quant

```

结果：

```
$ less -S quant.sf

Name    Length  EffectiveLength TPM     NumReads
ENST00000375213.8       3859    3576.267        16.992657       1048.590
ENST00000338983.7       7179    5823.482        10.274648       1032.439

```

从ENCODE下载的两个数据来看，在K562细胞系中，ZAKα的表达值比ZAKβ的要高一些。

  [1]: http://static.zybuluo.com/sherking/4u778jcsmb4phbvtj1n4uugo/image_1eg2om6bim3s2ude3dlk33eg9.png
  [2]: http://static.zybuluo.com/sherking/2qkw7hhgvzpxo4x6dgp9mc3d/image_1eg2omu413j3tf0c53163oo7hm.png
  [3]: http://static.zybuluo.com/sherking/lm306h0zlg3uzoa2c12ob8hn/image_1eg2qpqop2am1agc18h611c9if813.png