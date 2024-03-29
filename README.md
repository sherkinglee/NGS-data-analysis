﻿# **[NGS-data-analysis](https://sherkinglee.github.io/NGS-data-analysis/)**
## **背景**
2020年注定是不平凡的一年，寒假回家过年还是博士三年级，眼看快四年级了仍然躺在家中，面对网络的卡顿捉襟见肘。合作项目的实验室因为疫情的影响已经暂时关闭，MC的那篇文章估计还会拖很久；自己的文章被拒之后也一直没有接着往下投；后来想的两个项目因为都涉及到实验验证，自己所做的东西确实很有限。所以想趁着疫情恢复期，把自己之前一直想做但是没做的事情给做了，这也就是这个项目产生的背景了。
这个项目建立的初衷就是想把自己在博士期间所做过的数据分析简单整理一下，包括数据分析的流程、分析的脚本、绘图技巧、遇到的问题等等，这样方便自己以后再次遇到同样问题的时候不需要重复造轮子，直接拿过来用就好了。期间，可能不会涉及太多原理性的东西，毕竟只是自己用来查阅的一个资料，暂时不会对外开放，所以大部分还是会以脚本的形式出现，这样会方便我自己的复制粘贴，后面如果有机会的话会逐渐完善原理部分的介绍，尽可能形成一个系统的参考资料。

## **设计**
到目前为止，我接触的项目都是以二代测序数据分析为主，主攻的项目还是以核糖体图谱技术（ribosome profiling）为主的翻译调控。所以接触最多的就是Ribo-seq数据分析，期间也分析过RNA-seq、小RNA-seq、ChIP-seq、CLIP-seq以及单细胞转录组数据的分析，但是除了Ribo-seq数据分析之外，其他的都不是很深入。既然主要的分析都集中在NGS的数据上，那么我会先从NGS数据的下载，前期数据准备包括质控、去接头、去rRNA污染、mapping等开始，然后过渡到后续的差异表达分析、数据的绘图展示、功能富集分析、核酸和蛋白序列分析等等，最后再添加自己学习的一些其他类型的NGS数据分析内容。
## **目的**
目的很简单，就是不想再遇到同样问题的时候再重新设计脚本，直接拿过来用就好了，避免不必要的麻烦。

## **框架**

### **[前期数据的获取和质控到mapping](https://github.com/sherkinglee/NGS-data-analysis/blob/master/01.pre-processing/pre-processing.md)**

1） 上游数据的获取（GEO）数据的下载

2） 数据的质控（fastqc等）

3） adapter的去除

4） 去除rRNA的污染

5） [mapping （hisat2, STAR, bowtie, bowtie2 等不同mapping软件， 单端和双端的差别）](https://github.com/sherkinglee/NGS-data-analysis/blob/master/01.pre-processing/Mapping.md)


### **[拿到bam文件之后的数据分析](https://github.com/sherkinglee/NGS-data-analysis/blob/master/02.differential-expressioin-analysis/DEAnalysis.md)**

1）counts的计算，htseq-count， featureCounts等

2）标准化方式：RPKM, RPM, TPM，DESeq2标准化

3）数据的合并（R, bash, python)

4）差异分析 （DESeq2, edgeR, limma, xtail)

### **[拿到差异基因之后的数据分析，可视化和功能富集分析等](https://github.com/sherkinglee/NGS-data-analysis/blob/master/02.differential-expressioin-analysis/visualization.md)**

1） [火山图的绘制](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Volcano_plots.md)

2） [热图的绘制](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Heatmaps.md)

3） [GO和KEGG以及GSEA的数据分析](https://github.com/sherkinglee/NGS-data-analysis/tree/master/03.functional-analysis)

### **[其他](https://github.com/sherkinglee/NGS-data-analysis/blob/master/04.sequence-features-analysis)**

1）[motifs分析（MEME, Homer)](https://github.com/sherkinglee/NGS-data-analysis/blob/master/04.sequence-features-analysis/motifs.md)

2）[RNA二级结构： RNAfold, RNAstructure](https://github.com/sherkinglee/NGS-data-analysis/blob/master/04.sequence-features-analysis/RNASecondaryStructures.md)

3）[蛋白质二级结构预测](https://github.com/sherkinglee/NGS-data-analysis/blob/master/04.sequence-features-analysis/ProteinSecondaryStructures.md)

4）GC分析

5）UTR获取

6）保守性分析

7）[转录本定量---Salmon](https://github.com/sherkinglee/NGS-data-analysis/tree/master/04.sequence-features-analysis/isoformQuantification.md)

---

### **[常见图表绘制](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots)**


1）[火山图:volcano plot](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Volcano_plots.md)

2）[韦恩图:venn plot](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Venn_plots.md)

3）[Metagene plot](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Metagene_plots.md)

4）[累积密度分布图: cumulative distribution plot](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Cumulative_distribution_plots.md)

5）[热图:heatmap](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Heatmaps.md)

6）[核型图](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/KayotypePlot.md)

N）[将Rplots导出为PPT格式](https://github.com/sherkinglee/NGS-data-analysis/blob/master/05.plots/Rplots2PPT.md)

---

### **这些基本的流程弄完之后，再按照项目来整理**

1）RNA-seq数据分析流程

2）小RNA数据分析流程（miRNA数据分析）

3）[Ribo-seq数据分析流程](https://github.com/sherkinglee/NGS-data-analysis/blob/master/06.projects/Ribo-seq.md)

4）表观组数据分析流程（ATAC-seq，CHIP-seq,])

5）CLIP数据分析流程

6）单细胞转录组数据分析流程

7）单细胞空间转库组数据分析流程

8）miRNA的预测和targets的预测等

9）[MaxQuant分析pSILAC蛋白质组学数据](https://github.com/sherkinglee/NGS-data-analysis/tree/master/06.projects/MaxQuant-pSILAC.md)

---
### **[其他](https://github.com/sherkinglee/NGS-data-analysis/blob/master/07.others/)**

1）[snakemake流程搭建初步学习](https://github.com/sherkinglee/NGS-data-analysis/blob/master/07.others/snakemake.md)

+   [Ribo-seq数据分析流程](https://github.com/sherkinglee/NGS-data-analysis/blob/master/07.others/Ribo-seq-snakemake.py)
+ [RNA-seq数据分析流程](https://github.com/sherkinglee/NGS-data-analysis/blob/master/07.others/RNA-seq-snakemake.py)
+ [上述Ribo-seq和RNA-seq每步流程日志处理python脚本](https://github.com/sherkinglee/NGS-data-analysis/blob/master/07.others/summary.py)