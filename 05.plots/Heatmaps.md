## **Heatmaps**

热图是生物信息学领域常见的图表之一，在做二代测序数据分析差异分析之后往往需要进行热图等的绘制。热图有很多绘制方法，可以用R语言绘制，也可以用Python绘制，还可以用Excel进行绘制。在R语言里，用pheatmap，gplots, ggplot2, ComplexHeatmap等等绘制过热图，各有优缺点，选择一个就好；Python用seaborn和maplotlib也可以进行热图的绘制。Excel可以编写VBA程序绘制热图，也可以根据色阶进行绘制。

### **R: pheatmap**

#### **pheatmap安装**
```
install.packages(“pheatmap”) #安装pheatmap包
library(pheatmap) #加载pheatmap包
head(data_selected)
```
![image_1eeuq48d8es016ik17n41i3c64c9.png-193.1kB][1]

#### **pheatmap绘图**

```
library(pheatmap)
library(RColorBrewer)
library(circlize)
annotation_col=data.frame(type=factor(conditions$type))
rownames(annotation_col) <- conditions$conditions

myc <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
x<- pheatmap(data_selected,color = rev(myc),border_color = NA,scale = "none",cluster_rows=T,cluster_cols = T,clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",clustering_method = "ward.D2",show_rownames = F,display_numbers = F,annotation_col = annotation_col)
```

![image_1eeuq8o9s1ft1q611vd5vh515tdm.png-287.9kB][2]
  
  
#### **pheatmap获取聚类后的矩阵**
 
```
order_row = x$tree_row$order  #记录热图的行排序
order_col = x$tree_col$order    #记录热图的列排序
ordered_data = data.frame(data_selected[order_row,order_col])   # 按照热图的顺序，重新排原始数据
ordered_data = data.frame(rownames(ordered_data),ordered_data,check.names =F)  # 将行名加到表格数据中

```

--- 

### **R: gplots**

#### **gplots安装**

```
install.packages("gplots")
library(gplots)
?heatmap.2()
```

#### **gplots绘图**

```
library(gplots)
library(RColorBrewer)
myc <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
myDist=function(x){dist(x,method = "euclidean")}
myhclust = function(x){hclust(x,method="ward.D2")}
annotation_col=data.frame(type=factor(conditions$type))
rownames(annotation_col) <- conditions$conditions
head(annotation_col)
	
#type
#<fctr>
#FP-HS-2hrs-GFP(GSE32060)	Elongation			
#Hspa1a-HS-2hrs-Hspa1a(GSE32060)	Elongation			
#TSC2(+/-)-TSC2(+/+)(GSE78959)	Initiation			
#TSC2(-/-)-TSC2(+/+)(GSE78959)	Initiation			
#AZC-MG132-untreated(GSE71763)	Elongation			
#AZC-MG132-untreated-rRNA-del(GSE71763)	Elongation
colorsbar <- unlist(lapply(annotation_col$type,function(x){if(x=="Initiation")"#FF0000"else if (x=="Elongation")"#0000FF" else "#00FF00"}))

## expression change
y <- heatmap.2(as.matrix(data_selected),col=rev(myc),scale="none",trace="none",Rowv =TRUE,Colv=TRUE,dendrogram = "both",density.info="none",cexRow = 0.8,cexCol =0.8,margins = c(10,10) ,key.xlab="D(polarity)",distfun = myDist,hclustfun = myhclust,labRow = F,ColSideColors = as.vector(colorsbar))
legend(0.85,0.5, c("Initiation", "Elongation","RocA"), col=c("#FF0000", "#0000FF","#00FF00"), cex=0.8,pch = 19,border = F)

```

![image_1eev0ue241e1t1usc1lq418mq1k7413.png-177kB][3]

#### **gpots获取聚类后的数据矩阵**

```
ordered_data <- data_selected[rev(y$rowInd), y$colInd]
```


### **R:ComplexHeatmap**

#### **ComplexHeatmap绘图**
```
library(RColorBrewer)
library(circlize)
# library(export)
library(ComplexHeatmap)
annotation1 <- HeatmapAnnotation(df=data.frame(type=conditions$type),show_legend = T,col = list(type=c(Elongation="pink",Initiation="lightblue",RocA="lightgreen")),show_annotation_name = F,annotation_legend_param = list(title="Type",legend_direction="vertical"))

z<-Heatmap(as.matrix(data_selected),col = colorRamp2(c(-0.5,0,0.5), rev(colorRampPalette(brewer.pal(5, "RdBu"))(3))),cluster_rows = T,cluster_columns = T,clustering_distance_columns = 'euclidean',clustering_method_columns ="ward.D2",clustering_method_rows = "ward.D2",clustering_distance_rows = "euclidean",heatmap_legend_param = list(title="D(polarity)",legend_direction="vertical"),top_annotation = annotation1,column_names_gp = gpar(fontsize=8))

```

![image_1eev1q1dj147omqa17ddrf01q0e1g.png-264.5kB][4]

#### **complexHeatmap获取聚类之后的结果**

```
ordered_data = data_selected[row_order(z)[[1]],column_order(z)]
```

### **Python**

```
import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

fig=plt.figure()
sns.set(font_scale=0.2)
ax=sns.clustermap(data_selected,metric='euclidean',method='ward',row_cluster=True,col_cluster=True,
                vmin=-1,vmax=1,center=0,cmap = 'RdBu_r')
plt.subplots_adjust(left=0.08, bottom=0.205, right=0.835, top=0.93, wspace=0.2, hspace=0.2)
plt.show()
```
![image_1eev2fe5r1v3i47l135pm4mepa1t.png-318.8kB][5]


### **Excel**

```
Sub ChangeCellColor()

    For Each c In ActiveSheet.Range("C2:W21")
        If c.Value < 2 Then
            step = (c.Value - 0) / 2
            c.Interior.Color = RGB(217 + (254 - 217) * step, 217 + (251 - 217) * step, 217 + (212 - 217) * step)
        End If
        If (c.Value >= 2) And (c.Value < 3) Then
            step = (c.Value - 2) / 1
            c.Interior.Color = RGB(254 + (254 - 254) * step, 251 + (228 - 251) * step, 212 + (146 - 212) * step)
        End If
        If (c.Value >= 3) And (c.Value < 4) Then
            step = (c.Value - 3) / 1
            c.Interior.Color = RGB(254 + (254 - 254) * step, 228 + (196 - 228) * step, 146 + (76 - 146) * step)
        End If
        If (c.Value >= 4) And (c.Value < 6) Then
            step = (c.Value - 4) / 2
            c.Interior.Color = RGB(254 + (254 - 254) * step, 196 + (154 - 196) * step, 76 + (42 - 76) * step)
        End If
        If (c.Value >= 6) And (c.Value < 10) Then
            step = (c.Value - 6) / 4
            c.Interior.Color = RGB(254 + (217 - 254) * step, 154 + (95 - 154) * step, 42 + (14 - 42) * step)
        End If
        If (c.Value >= 10) And (c.Value < 20) Then
            step = (c.Value - 10) / 10
            c.Interior.Color = RGB(217 + (153 - 217) * step, 95 + (52 - 95) * step, 14 + (4 - 14) * step)
        End If
        If c.Value >= 20 Then
            c.Interior.Color = RGB(163, 55, 4)
        End If
        'c.Interior.Color = RGB(255, 255, 255)
    Next c
End Sub


```

![image_1eev3pice1u4d1l541ga61c811k072a.png-73.9kB][6]


  [1]: http://static.zybuluo.com/sherking/556isf74b45mspdao2gxve45/image_1eeuq48d8es016ik17n41i3c64c9.png
  [2]: http://static.zybuluo.com/sherking/651zgoywfm87px2qteaa8x55/image_1eeuq8o9s1ft1q611vd5vh515tdm.png
  [3]: http://static.zybuluo.com/sherking/g6pzrawzajp15252c4qc2spr/image_1eev0ue241e1t1usc1lq418mq1k7413.png
  [4]: http://static.zybuluo.com/sherking/7mf2a0cfek84anj1aai24dji/image_1eev1q1dj147omqa17ddrf01q0e1g.png
  [5]: http://static.zybuluo.com/sherking/knhso7l0upl0e3hpcxw1gnzy/image_1eev2fe5r1v3i47l135pm4mepa1t.png
  [6]: http://static.zybuluo.com/sherking/h6jpljundyajexh3w08zcpq8/image_1eev3pice1u4d1l541ga61c811k072a.png