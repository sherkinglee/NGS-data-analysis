## **基因的差异表达分析之后的可视化**
这里的可视化分析不是一般常用绘图脚本，主要还是承接上一个文档的差异表达分析之后需要做的一些常规绘图，主要包括火山图的绘制，热图的绘制等等。至于常见图标的绘制，我们会有专门的一个主题来进行介绍。

## **火山图绘制**
火山图的绘制有很多种方式，这里我们选用ggplot2进行绘制

```
## 默认输入起始文件是做完差异表达之后的文件， 构建画图用的数据
log2FC <- xtail_results_dataframe$log2FC_TE_final
padjust <- xtail_results_dataframe$pvalue.adjust
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")
x_lim <- max(log2FC,-log2FC)

## 载入ggplot2等package
library(ggplot2)
library(RColorBrewer)

## 绘图
p <- ggplot(data_for_plot,aes(log2FC,-1*log10(padjust),color = sig))+geom_point()+xlim(-5,5)+labs(x="log2TE(Torin/Ctrl)",y="-log10(p.adjust)")

p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,15)
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
```
![image_1e6tval2a1kiq7467rj1va2bbg13.png-108.7kB][1]


## **热图绘制**
热图的绘制也有很多方式，我们这里选用gplots以及ComplexHeatmap

第一种： gplots绘图

```
## 安装gplots
install.packages('gplots')

## 载入gplot
suppressMessages(library(gplots))
suppressMessages(require(RColorBrewer))
myc <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

## expression change
heatmap.2(as.matrix(all_TOP),col=rev(myc),scale="none",trace="none",Rowv =TRUE,Colv=TRUE,density.info="none",cexCol =0.8,margins = c(10,10) ,key.xlab="log2FC")

## 其它参数详见
?plots()
```
![image_1e6u018mpt7213nfaq8s2316lh1t.png-104.6kB][2]

如果需要添加分类bar则可以这样

```
colors <- unlist(lapply(annotation_col$condition,function(x){if(x=="control")"#FF0000"else"#0000FF"}))
heatmap.2(as.matrix(DE_genes_counts_normalized),scale = "column",trace = "none",col = greenred(75),density.info = "none",cexCol = 1,labRow = FALSE,dendrogram = "column",ColSideColors = as.vector(colors))
```
![image_1e6vbjaor13he1n59t19n2ivcv9.png-35.2kB][3]

第二种： ComplexHeatmap

```
## 载入工具包
library(RColorBrewer)
library(circlize)
library(export) ## 用于输出ppt格式的图表
library(ComplexHeatmap)

## 绘图
z<-Heatmap(as.matrix(up_genes_MS_RPF_com[,c(1,2)]),col = colorRamp2(c(-2,0,2), rev(colorRampPalette(brewer.pal(5, "RdBu"))(3))),cluster_rows = T,cluster_columns = F,clustering_distance_columns = 'euclidean',clustering_method_columns ="ward.D2",heatmap_legend_param = list(title="log2FC",legend_direction="vertical"))
z
```
![image_1e6u0fbbt1iosj0nuhh1uvs15ii9.png-71.3kB][4]


第三种： pheatmap

```
# plot heatmap using pheatmap
suppressMessages(library(pheatmap))
DE_genes_counts_normalized <- log2(DE_genes_counts_normalized)
annotation_col <- as.data.frame(colData(dds)[,c("sample","condition")])
pheatmap(DE_genes_counts_normalized,show_rownames = FALSE,cutree_rows = 4,cutree_cols = 2,annotation_col = annotation_col)
```
![image_1e6vbmm26c9s1jn414fg1l3o12h4m.png-73.9kB][5]


  [1]: http://static.zybuluo.com/sherking/m13toiu8pat9bbaxhjy7bg0k/image_1e6tval2a1kiq7467rj1va2bbg13.png
  [2]: http://static.zybuluo.com/sherking/yu7wseizfo1t1lzebu3zxj6c/image_1e6u018mpt7213nfaq8s2316lh1t.png
  [3]: http://static.zybuluo.com/sherking/09nvljsp038b6w793mf54cgm/image_1e6vbjaor13he1n59t19n2ivcv9.png
  [4]: http://static.zybuluo.com/sherking/0g4wtqfejd3km5fbqvvb8wsn/image_1e6u0fbbt1iosj0nuhh1uvs15ii9.png
  [5]: http://static.zybuluo.com/sherking/eersbokwtacaqsq2522zdslu/image_1e6vbmm26c9s1jn414fg1l3o12h4m.png