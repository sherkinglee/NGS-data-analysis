## **差异表达分析**



序列比对完之后需要进行read counts的计算，然后进行差异表达分析以及后续的可视化处理和功能分析等等。我刚开始接触到差异分析表达是2017年在清华大学鲁志老师实验室学习到的，当时还是生物信息小白一个，先做的是芯片数据的分析处理，当时用到的工具包就是limma，后面进一步的轮转学习开始逐步接触更多的数据分析工具，其中RNA-seq数据用的最多的就是DESeq2，edgeR等R包，目前Ribo-seq数据分析的也常用我们实验室开发的xtail进行翻译效率（TE）的差异分析，当然也可以继续使用DESeq2进行Ribosome footprints(RPF)的差异分析。这部分我主要从reads counts的计算， 数据的合并，标准化方式的介绍以及不同packages的使用进行介绍。

---
## **测序读段的计数（read counts的计算）**
目前计数的工具很多，根据是否需要mapping可以分为alignment-based和alignment-free两种，具体的可以网络上的一篇[博客](https://www.plob.org/article/11504.html)。表达数据的定量，分为三个层次： gene level, transcript level and exon level. 不同的level也有不同的计算方式，这里我们以经典的HTseq对基因level的表达进行定量。

第一步： [HTseq](https://htseq.readthedocs.io/en/master/)下载和安装

```
## python 下载
pip install HTSeq
```

第二步： htseq-count计数

```
#!/bin/bash
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE70211/11.read_counts
Ref=/workdata/LabMember2/lifj/lifj/data/Reference/human
BamDir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE70211/07.STAR

for i in SRR20759{30..34};do

htseq-count -f bam -s yes -i gene_name -t exon -m intersection-strict $BamDir/${i}_STAR/${i}.Aligned.sortedByCoord.out.bam $Ref/Homo_sapiens.GRCh38.88.gtf > $workdir/${i}_intersection-strict_exon_STAR.counts

具体参数设置详看htseq-count -h
```
第三步：查看htseq-count的输出结果

```
## htseq-count的结果大致如下，第一列为ID，可以是gene_name，可以是gene_id等等，第二列是counts数
5S_rRNA 1
5_8S_rRNA       0
7SK     0
A1BG    15
A1BG-AS1        39
A1CF    4
A2M     4
A2M-AS1 0
A2ML1   0
A2ML1-AS1       5
A2ML1-AS2       0
A2MP1   0

```
## **不同样本数据的合并（生成counts 矩阵的过程）**

如上述，拿到htseq-count的计数结果之后如果需要进行后续分析，最好需要把不同样本的counts文件合并成一个矩阵这样方便统一标准化和处理。在我的分析经验中有很多种不同的方式可以把表达文件合并，bash/R/python都可以。

第一种： bash下合并表达文件: paste+awk。

```
## 利用paste和awk计算
paste *.counts | awk '{printf $1 "\t";for(i=2;i<=NF;i+=2) printf $i"\t";printf "\n"$i}' > test/bash_counts.txt

## 结果如下：缺点就是没有列名，需要手动添加，或者用sed添加一行列名： sed '1igene_name\tsample1\tsample2\tsample3\tsample4\t'
5S_rRNA 1       0       0       0       1
5_8S_rRNA       0       0       0       0       0
7SK     0       0       0       0       0
A1BG    12      12      14      16      15
A1BG-AS1        46      37      50      42      39
A1CF    3       8       1       3       4
A2M     5       3       2       3       4
A2M-AS1 0       0       0       0       0
A2ML1   1       0       0       0       0
A2ML1-AS1       7       4       5       7       5
A2ML1-AS2       0       0       0       0       0
A2MP1   0       2       0       0       0
```

第二种： R语言合并文件

```
## 第一步：一次读入每个表达数据文件
sample1 <- read.table("sample1.txt",stringsAsFactors = F,header = T,check.names = F,sep="\t")
sample2 <- read.table("sample2.txt",stringsAsFactors = F,header = T,check.names = F,sep="\t")
...
sample5 <-  read.table("sample5.txt",stringsAsFactors = F,header = T,check.names = F,sep="\t")

## 第二步： 编写合并函数
## when the element in the list only has two columns, one is name, and the other is the value, the function is useful
multimerge <- function(data=list()){
  if(length(data)< 1){
    stop("Input error")
  }
  else if(length(data)==1){
    tem=data[[1]]
    return(as.data.frame(tem))
  }
  else{
    tem=data[[1]]
    tem=tem[,-2]
    for(i in data){
      tem=data.frame(tem,i[,2])
    }
    rownames(tem) <- as.character(tem[,1])
    tem=tem[,-1]
  }
  return(tem)
}

## 第三步： merge 数据
data <- multimerge(list(sample1,sample2,sample3,sample4,sample5))
colnames(data) <- c('sample1','sample2','sample3','sample5')

```

第三种： python合并数据

```
略
```

## **标准化方式（RPKM，RPM, TPM, DESeq2等）**
在NGS数据处理过程中需要考虑很多能够影响不同样品之间表达值比较的因素，包括基因长度、测序深度等等。目前有很多标准化的方式可以在某种程度上消除这种因素的影响，包括RPM, RPKM, TPM, DESeq2标准化方式等等。

第一种： RPM/CPM (Reads/Counts of exon model per Million mapped reads)

```
## 计算公式
RPM/CPM= 10^6 * (read counts of a gene)/(total mapped reads of all genes)

## R 语言中实现
data_RPM <- apply(data, 2, function(x){10^6*(x/sum(x))})
```

第二种： RPKM （Reads Per Kilobase of exon modelper Million mapped reads）（单端RPKM, 双端：FPKM)

在这里最主要的参数就是测序深度和基因长度的定义了，我一般定义为测序深度就用有效mapped的reads总数，即htseq-count计算出来的所有基因上的reads总数，而基因的长度可以有两种使用策略，一种是直接用每个基因最长的转录本长度进行计算，另一个策略就是把每个基因的exon提取出来，然后去掉overlap的部分长度，最后求和得到就是这个基因所有exon长度的总和。第二种我其实用的还不多，主要还是用第一种，用每个基因最长转录本长度来计算，具体的计算方式可以用python计算即可，详情可以参见[RiboMiner](https://github.com/xryanglab/RiboMiner)中的[OutputTranscriptInfo](https://github.com/xryanglab/RiboMiner/blob/master/RiboMiner/OutputTranscriptInfo.py)函数。
```
## 计算方式
RPKM=10^9*(read counts of a gene / (total mapped reads of all genes * length of a gene))

## R 语言中实现 (转录本最长长度)
data <- multimerge(list(sample1,sample2,sample3,sample4,sample5))
colnames(data) <- c('sample1','sample2','sample3','sample5')
trans_length <- read.table("longest.trasn.info.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
## use common gene of the longest and the htseq results
略
## RPKM 计算
total_counts <- colSums(data)
data_RPKM <- t(do.call(rbind,lapply(1:length(total_counts),function(x){
10^9*data[,i]/(trans_length$CDS_length*total_count)
})))

## R 语言实现 （GemomicFeatures) （参见https://www.bioinfo-scrounger.com/archives/342/）
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("hg38.gtf",format="gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- lapply(exons_gene,function(x){sum(width(reduce(x)))})
## 得到exon的总长度之后就可以按照上面的公式进行计算RPKM值了

```

第三种： TPM（Transcripts Per Kilobase of exonmodel per Million mapped reads）,推荐用RSEM计算得到。此处参考[博文](www.bio-info-trainee.com/2017.html)。


![image_1e6qqf21h1930nic1g2flf2bji9.png-123.8kB][1]

```
## 按照上述方法计算得到RPKM值之后，按照这个图片的公式计算TPM值

data_TPM <- apply(data_RPKM,2,function(x){10^6*x/sum(x)})
```

第四种： DESeq2标准化方式，参见[DESeq2官网](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)。部分参见[微信公众号](https://mp.weixin.qq.com/s/Vmhx_TGxNkQzkekf93Xl4w)。

```
##安装DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## DESeq2标准化
library(DESeq2)
normalize_factor <- estimateSizeFactorsForMatrix(data)
data_normed <- t(t(data)/normalize_factor)

## 或者先进行DESeq2的标准化，然后输出标准化的数据
normalized_counts <- counts(dds, normalized=TRUE)
```

## **差异表达分析（limma, edgeR, DESeq2, xtail等）**

### **[limma使用](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)**

```
## 安装limma
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("limma")

## 载入limma
library(limma)

## 读入表达矩阵，一般会将低表达的数据过滤掉，这里都省略
data <- read.table("data.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(data) <- data[,1]
data <- data[,-1]

## 设计矩阵
groups <- c("Normal","Normal","Tumor","Tumor")
groups <- factor(groups,levels = c("Normal","Tumor"))
designs <- model.matrix(~0+groups)
colnames(designs) <- c("Normal","Tumor")
rownames(designs) <- c("Normal-1","Normal-2","Tumor-1","Tumor-2")

## 数据拟合
fit <- lmFit(data,designs)

## 构建对照矩阵,设定谁比谁
contrast.matrix<-makeContrasts(Tumor-Normal,levels = designs)
contrast.matrix
#        Contrasts
#Levels   Tumor - Normal
#  Normal             -1
#  Tumor               1

fit2 <- contrasts.fit(fit,contrast.matrix)

## 进行经验贝叶斯检验
fit2 <- eBayes(fit2)

## 差异表达结果
diff_results <- topTable(fit2, coef=1, n=Inf,adjust.method ='fdr')

```




### **[edgeR使用](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html)**

```
## 安装edgeR,R>
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")

## 读入表达矩阵
library(edgeR)
data <- read.table("data.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(data) <- data[,1]
data <- data[,-1]

## 构建DGElist对象
groups <- factor(c("Normal","Normal","Tumor","Tumor"))
y <- DGEList(counts=data,group=groups)
y

## 过滤掉低表达数据
keep<-rowSums(cpm(y)>1)>=2
y<-y[keep,,keep.lib.sizes=FALSE]
y

## 标准化TMM
y<-calcNormFactors(y)
y

## 构建实验设计矩阵
designs <- model.matrix(~0+groups)

## 评估离散度
y<-estimateDisp(y,designs,robust=TRUE) 
plotBCV(y)

## 差异分析，可以采用不同的拟合策略
fit<-glmQLFit(y,designs,robust=TRUE) # 似然比F检验
qlf<-glmQLFTest(fit)

## 也可以
fit <- glmFit(y, designs) #根据分布进行数据处理
lrt <- glmLRT(fit, coef=2) ＃返回差异表达的基因，计算统计量

## 具体问题实际操作过程再看

```


### **[DESeq2使用](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html)**

```
## DESeq2安装
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

## DESeq2使用
library(DESeq2)

## 读入表达数据
data <- read.table("data.txt",sep="\t",header=T,stringsAsFactors=F,check.names=F)
rownames(data) <- data[,1]
data <- data[,-1]

## 构建信息矩阵
colData <- data.frame(samples=colnames(data),condition=rep(c("Tumor","Normal"),each=2))
rownames(colData) <- colData$samples

## 转换数据格式
data_dds <- DESeqDataSetFromMatrix(countData = data,colData = colData,design=~condition)

## 过滤掉低表达数据
data_dds <- data_dds[rowSums(counts(data_dds)) >=10]

## 构建分组因子
data_dds$condition <- factor(data_dds$condition,levels = c("Tumor","Normal"))

## 差异分析
data_dds <- DESeq(data_dds)
data_res <- results(data_dds)

## 去掉NA值
data_res <- na.omit(data_res)

## 获取标准化后的表达数据
normalized_counts <- counts(data_dds, normalized=TRUE)

```



### **[Xtail使用](https://github.com/xryanglab/xtail/blob/master/vignettes/)**

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

  [1]: http://static.zybuluo.com/sherking/vmhx330fwwi00cpxk9nz75x8/image_1e6qqf21h1930nic1g2flf2bji9.png