## 染色体核型分析 ##

 

染色体核型分析也是生物信息中常见的一个分析，绘图的软件包也有很多，这里主要介绍两种**[RIdeogram](https://mran.microsoft.com/web/packages/RIdeogram/vignettes/RIdeogram.html)**和**[karyoploteR](https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotRegions/PlotRegions.html)**

### **[RIdeogram](https://mran.microsoft.com/web/packages/RIdeogram/vignettes/RIdeogram.html)**

+ 下载和安装

```
install.packages("RIdeogram")
```

+ 使用

```
library(RIdeogram)
mouse_karyotype <- read.table("mouse_karyotype.txt",sep="\t",stringsAsFactors = F,header = T,check.names = F)
gene_density <- GFFex(input = "Mus_musculus.GRCm38.87.gtf", karyotype = "mouse_karyotype.txt", feature = "gene", window = 1000000)
head(mouse_karyotype)
  Chr Start       End
1   1     0 195471971
2   2     0 182113224
3   3     0 160039680
4   4     0 156508116
5   5     0 151834684
6   6     0 149736546

head(gene_density)

  Chr   Start   End Value
1   1       1 1e+06     0
2   1 1000001 2e+06     0
3   1 2000001 3e+06     0
4   1 3000001 4e+06    17
5   1 4000001 5e+06    25
6   1 5000001 6e+06    10


## prepare lncRNA density
lncRNA_density <- data.frame(Type=zz_transID_clusters$clusters,Shape=rep(c("circle","box","triangle"),times=c(nrow(polyA_raw_counts_common_cluster1),nrow(polyA_raw_counts_common_cluster2),nrow(polyA_raw_counts_common_cluster3))),Chr=transInfo[match(zz_transID_clusters$transID,transInfo$transcriptID),2],Start=transInfo[match(zz_transID_clusters$transID,transInfo$transcriptID),3],End=transInfo[match(zz_transID_clusters$transID,transInfo$transcriptID),4],color=rep(c("6a3d9a","33a02c","ff7f00"),times=c(nrow(polyA_raw_counts_common_cluster1),nrow(polyA_raw_counts_common_cluster2),nrow(polyA_raw_counts_common_cluster3))))


ideogram(karyotype = mouse_karyotype,output = "mouse_chromosome2.svg",overlaid = gene_density,label = lncRNA_density,label_type = "marker")
convertSVG("mouse_chromosome2.svg",file = "mouse_chromosome2", device = "png")
```

![image_1fcd5uuee1pba1u54d03rth1e6a9.png-258.4kB][1]


### **[karyoploteR](https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotRegions/PlotRegions.html)**

+ 下载和安装

```
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("karyoploteR")
```

+ 使用

```
library(karyoploteR)
library(rtracklayer)

gtf.file <- "lncRNA_3624_common2.gtf"

features <- import(gtf.file)
table(features$type)
genes <- features[features$type=="transcript"]


pp <- getDefaultPlotParams(plot.type=2)
pp$data1outmargin <- 100
pp$data2outmargin <- 100
pp$topmargin <- 450
kp <- plotKaryotype(genome="mm10", ideogram.plotter = NULL,plot.type=2, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddMainTitle(kp, "mm10 with novel lncRNAs", cex=1.1)
kpPlotRegions(kp, data=genes[strand(genes)=="+"], avoid.overlapping = FALSE, col="deepskyblue",lwd=1.5)
kpPlotRegions(kp, data=genes[strand(genes)=="-"], avoid.overlapping = FALSE, col="gold", data.panel=2,lwd=1.5)
kpAddLabels(kp, "(+)", cex=0.8, col="#888888")
kpAddLabels(kp, "(-)", data.panel=2, cex=0.8, col="#888888")
```

![image_1fcd61se1aa41i3u1qvb7kg1gm8m.png-113kB][2]


  [1]: http://static.zybuluo.com/sherking/5f9uky1vp5ihoe173my8aa0d/image_1fcd5uuee1pba1u54d03rth1e6a9.png
  [2]: http://static.zybuluo.com/sherking/3lw4ucu2zcks77lkoonoc28h/image_1fcd61se1aa41i3u1qvb7kg1gm8m.png