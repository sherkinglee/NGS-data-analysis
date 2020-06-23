## **Volcano plot**



火山图的绘制一般用于差异基因表达之后用于展示在某种条件下上下调的基因情况。关于火山图，几乎每个生信人员都会绘制，这里我一般用ggplot2进行绘制，当然也可以结合ggpubr, ggtheme等R package进行绘制。

### **火山图： ggplot2**

+ 数据准备
进行火山图绘制所需要的数据一般都是进行差异表达之后的数据，即DESeq2,limma, xtail等差异表达的输出文件，包括log2FC, pvalue, qvalue等值。大概数据类型如下：

```
         baseMean log2FoldChange     lfcSE      stat       pvalue         padj
EIF3E    526.7062      -3.176044 0.1544955 -20.55752 6.592013e-94 8.481284e-90
IFI6    7896.6262       2.312131 0.1138886  20.30169 1.242253e-91 7.991411e-88
BST2     927.4316       2.188719 0.1234603  17.72812 2.543834e-70 1.090966e-66
KLK5    1588.1299      -2.200539 0.1247072 -17.64565 1.099104e-69 3.535266e-66
KLK7     529.2680      -2.447886 0.1482081 -16.51655 2.789081e-61 7.176864e-58
SLC39A8  657.2581       1.927590 0.1237401  15.57773 1.031606e-54 2.212106e-51
```

+ R 包安装

```
install.packages("ggplot2")
```
+ 绘制火山图

```
library(ggplot2)
DE_results <- as.data.frame(res)
log2FC <- DE_results$log2FoldChange
padjust <- DE_results$padj

## prepare dataset for plot
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")

## volcano plot
x_lim <- max(log2FC,-log2FC)
theme_set(theme_bw())
p <- ggplot(data_for_plot,aes(log2FC,-1*log10(padjust),color = sig))+geom_point()+xlim(-5,5)+labs(x="log2TE(Treat/Ctrl)",y="-log10(p.adjust)")

p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))+ylim(0,15)
# p <- p  +guides(colour = FALSE) ## 选择是否需要图例
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text = element_text(size=15),legend.title = element_text(size=15))
p
```

![image_1ebg69q5ds2oo0u1phq1oqsnqb9.png-67.9kB][1]

+ 输出plot为ppt格式

```
## 利用export package
install.pacakge('export')
library(export)
graph2ppt(p,"volcano_plot.ppt")

## 但是export被下架，所以可以选择另外一个eoffice,结合ggplotify可以把任何R plot输出为PPT格式
install.package('eoffice')
library(eoffice)
topptx(g,'volvano_plot2.pptx')

## 如果是基础包画出的图形，可以用ggplotify转化为ggplot类型图标
g=as.ggplot(plot(x,y))
topptx(g,'scatter_plot.pptx')

```

+ 绘制火山图+带gene label

```
library(ggplot2)
library(ggrepel)
DE_results <- as.data.frame(res)
log2FC <- DE_results$log2FoldChange
padjust <- DE_results$padj

## prepare dataset for plot
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")
rownames(data_for_plot) <- rownames(DE_results)[which(DE_results$log2FoldChange%in%data_for_plot$log2FC)]

## 挑选上下调变化最显著的前10个基因，p排序
data_for_plot$label=NA
data_for_plot <- data_for_plot[order(data_for_plot$padjust,decreasing = F),]
data_for_plot$label[head(which(data_for_plot$sig=="up"),10)]<- rownames(data_for_plot)[which(data_for_plot$sig=="up")]
data_for_plot$label[head(which(data_for_plot$sig=="down"),10)]<- rownames(data_for_plot)[which(data_for_plot$sig=="down")]



## volcano plot
x_lim <- max(log2FC,-log2FC)
theme_set(theme_bw())
p <- ggplot(data_for_plot,aes(log2FC,-1*log10(padjust),color = sig))+geom_point()+xlim(-5,5)+labs(x="log2TE(Treat/Ctrl)",y="-log10(p.adjust)")

p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)  + geom_text_repel(aes(label = label), box.padding = unit(0.3, "lines"), point.padding = unit(0.4, "lines"), show.legend = F, size = 3)
p <- p +theme(panel.grid =element_blank())+theme(axis.line = element_line(size=0))
# p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text = element_text(size=15),legend.title = element_text(size=15))
p
```

![image_1ebgb72r21hc0h261c391n2j1l4vp.png-61.4kB][2]

### **火山图： ggpuber+ggthemes**

和直接用ggplot产生的图差不多

```
# install.packages('ggpubr')
# install.packages('ggthemes')
DE_results <- as.data.frame(res)
log2FC <- DE_results$log2FoldChange
padjust <- DE_results$padj

## prepare dataset for plot
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")
rownames(data_for_plot) <- rownames(DE_results)[which(DE_results$log2FoldChange%in%data_for_plot$log2FC)]

## 挑选上下调变化最显著的前10个基因，p排序
data_for_plot$label=NA
data_for_plot <- data_for_plot[order(data_for_plot$padjust,decreasing = F),]
data_for_plot$label[head(which(data_for_plot$sig=="up"),10)]<- rownames(data_for_plot)[which(data_for_plot$sig=="up")]
data_for_plot$label[head(which(data_for_plot$sig=="down"),10)]<- rownames(data_for_plot)[which(data_for_plot$sig=="down")]
data_for_plot$logP=-1*log10(data_for_plot$padjust)


## volcano plot
x_lim <- max(log2FC,-log2FC)

p <- ggscatter(data_for_plot,x="log2FC",y="logP",color="sig",palette = c("#0072B5","grey","#BC3C28"),size = 2,label = data_for_plot$label,repel = T,xlab = "log2TE(Treat/Ctrl)",ylab = "-log10(p.adjust)",font.label = 8) + theme_bw()+ geom_hline(yintercept=-log10(0.05),linetype=4)+ geom_vline(xintercept=c(-1,1),linetype=4) 


p <- p +theme(panel.grid =element_blank())+theme(axis.line = element_line(size=0))
# p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text = element_text(size=15),legend.title = element_text(size=15))
p

```

![image_1ebgcueaudt814s1cin1cok16t716.png-132.3kB][3]


### **交互式火山图： ggplot2+plotly**

```
library(plotly)
library(ggplot2)
DE_results <- as.data.frame(res)
log2FC <- DE_results$log2FoldChange
padjust <- DE_results$padj

## prepare dataset for plot
data_for_plot <- data.frame(log2FC=log2FC,padjust=padjust)
data_for_plot <- na.omit(data_for_plot)
data_for_plot$sig[(data_for_plot$padjust>0.05|data_for_plot$padjust=="NA")|(data_for_plot$log2FC>-1|data_for_plot$log2FC<1)] <- "NO"
data_for_plot$sig[(data_for_plot$padjust<=0.05&data_for_plot$log2FC>=1)] <- "up"
data_for_plot$sig[(data_for_plot$padjust<=0.05 & data_for_plot$log2FC<=-1)] <- "down"
summary(data_for_plot$sig=="up")
summary(data_for_plot$sig=="down")
summary(data_for_plot$sig=="NO")
rownames(data_for_plot) <- rownames(DE_results)[which(DE_results$log2FoldChange%in%data_for_plot$log2FC)]

data_for_plot$label<- rownames(data_for_plot)

## volcano plot
x_lim <- max(log2FC,-log2FC)
theme_set(theme_bw())
p <- ggplot(data_for_plot,aes(log2FC,-1*log10(padjust),color = sig))+geom_point(aes(text = label), alpha=0.8, size = 1)+xlim(-5,5)+labs(x="log2TE(Treat/Ctrl)",y="-log10(p.adjust)") 

p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4) 
p <- p +theme(panel.grid =element_blank())+
    theme(axis.line = element_line(size=0))
# p <- p  +guides(colour = FALSE) ## 选择是否需要图例
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20),legend.text = element_text(size=15),legend.title = element_text(size=15))
p <- ggplotly(p)
p
```

![image_1ebge3bl11l6n170idia1elr1ctm20.png-138.7kB][4]


  [1]: http://static.zybuluo.com/sherking/w321270dxi6w76qvju0xqkkf/image_1ebg69q5ds2oo0u1phq1oqsnqb9.png
  [2]: http://static.zybuluo.com/sherking/ula2izx8lkb50cmp9i1u6po8/image_1ebgb72r21hc0h261c391n2j1l4vp.png
  [3]: http://static.zybuluo.com/sherking/gpev64rlucr51bv87rfswdg5/image_1ebgcueaudt814s1cin1cok16t716.png
  [4]: http://static.zybuluo.com/sherking/c8t5suxbe589b8pe00cu38zy/image_1ebge3bl11l6n170idia1elr1ctm20.png