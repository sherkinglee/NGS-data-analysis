## **Motifs analysis**

常用的motifs分析软件主要包括[MEME](http://meme-suite.org/)和[Homer](http://homer.ucsd.edu/homer/)。[Homer](http://homer.ucsd.edu/homer/)被设计出来主要是分析ChIP-seq数据的，后面也可以用于其他peak file的处理和分析,包括GRO-Seq, RNA-Seq, DNase-Seq, Hi-C。我的分析经验有限，主要是用MEME和homer进行de novo motifs的identification. 其他经验可以参考网上的资料，其中包括清华大学[鲁志老师实验室的学习资料](https://lulab2.gitbook.io/teaching/part-iii.-ngs-data-analyses/5.motif/sequence_motif)

### **[Homer](http://homer.ucsd.edu/homer/)**

+ 软件下载详见[homer download](http://homer.ucsd.edu/homer/introduction/install.html).
+ Motif finding using fasta file

```
 findMotifs.pl targets.fa fasta OutputDirectory -fasta background.fa
 
 ## for example
 findMotifs.pl 374_up_CDS.fa human up_CDS -fasta ~/Reference/human/longest_CDS.fa

```

+ Motif finding using peak file

```
findMotifsGenome.pl peaks.txt ALF.fasta OutputResults/

## for example
findMotifsGenome.pl peaks.txt hg38 peakAnalysis -size 200 -len 8
```

+ peak annotation

```
## for example
peakFile=/Share2/home/lifj/Projects/01.RiboPipe/human/GSE70211/GSE85155_CLIP/00.peaks/GSE85155_G2_heIF4A1_hg19.narrowPeak
annotationFile=/Share2/home/lifj/Reference/human/hg19/Homo_sapiens.GRCh37.87.gtf

annotatePeaks.pl $peakFile hg19 -gtf $annotationFile  > annotationPeaks.txt

```


### **[MEME](http://meme-suite.org/)**

MEME-suit包括在线工具和单机版的本地工具，两种都可以用，本地工具的参数设置可以更加灵活一些，之前做过测试，在线工具和本地工具call出来的motifs有时候还是存在很大差别的，建议用Homer和MEME都用一下，看看那个更加合适一些。

+ 软件下载见[MEME Download](http://meme-suite.org/doc/download.html)
+ 第一步： 获取fasta序列。
    如果拿到的是基因组坐标，那么可以用bedtools getfasta获取对应的fasta序列，然后用MEME进行motifs 的identifacation。

+ 第二步： Motif finding using fasta file

```
## for example

meme $workdir/374_up_CDS.fa -dna -oc $workdir/up_CDS -nostatus -mod anr -nmotifs 10 -evt 1e+005 -minsites 2 -maxsites 500 -maxsize 10000000 -minw 5 -maxw 20 -revcomp
```


### **注意**

+ homer输出含有html的报告文件，包括de novo的motifs，也包括known motifs，如果觉得homer输出的logo不好看，可以提取其pwm文件，然后用R 包seqLogo进行绘制

```
## pwm 文件：
A	C	G	T
0.894	0.023	0.001	0.082
0.001	0.040	0.859	0.100
0.751	0.247	0.001	0.001
0.668	0.107	0.224	0.001
0.001	0.001	0.809	0.189
0.134	0.816	0.001	0.049
0.360	0.001	0.638	0.001
0.697	0.022	0.001	0.280

## seqLogo
library(seqLogo)
down_motif7 <- read.table("down_5UTR_motifs_7.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
down_motif7 <- t(down_motif7)
# pdf("down_5UTR_motif7.pdf")
seqLogo(down_motif7,ic.scale = F)
# dev.off()
```

+ MEME的输出结果里，理论上也会有html的报告文件，但是由于系统原因这个html文件可能不会生成，最后产生的只有不同motif的eps文件，如果需要可以利用PS或者AI转化成png或者pdf文件进行查看。
+ 其他关于motifs的分析也可以参见清华大学[鲁志老师实验室的学习资料](https://lulab2.gitbook.io/teaching/part-iii.-ngs-data-analyses/5.motif/sequence_motif)