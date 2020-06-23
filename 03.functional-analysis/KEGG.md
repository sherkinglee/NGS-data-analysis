## **KEGG pathway analysis**

KEGG分析在生物信息分析过程中也非常常见，目前有很多方法可以进行KEGG分析，包括使用[DAVID](https://david.ncifcrf.gov/), [Metascape](https://metascape.org/gp/index.html#/main/step1)以及[ClusterProfiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)等等。这里我们分别进行介绍。

### **[DAVID基本使用](https://david.ncifcrf.gov/)**

第一步：点击[Functional Annotation](https://david.ncifcrf.gov/summary.jsp)进入功能注释界面

![image_1e6v7iaqo1qg87lg14021etq1svn9.png-231.9kB][1]

第二步： 上传gene list
![image_1e6v7r2bc1ohuks1pdb1jn011l6m.png-122kB][2]

  
第三步： 选择gene list的ID类型，比如说gene_name(gene symbol), gene id ,transcript id等等

![image_1e6v7tsbl1eeetk7kcd1vsv1b9e13.png-160.3kB][3]


第四步： 选择gene list的类型，包括目标List和background list，最后提交上去就可以了。这里需要注意的是如果不提供background则默认以所有基因为背景，可以选择自己提供的gene list为背景。


![image_1e6v829pu35egbk1mqk13k7p2p1g.png-84kB][4]


第五步： 选择gene list的物种信息和功能注释的类别

![image_1e6vdr6bt1ikv1urf1oet1kgs039.png-118kB][5]


第六步： 查看结果

![image_1e6vs1glf11ii15gocq41i422331t.png-60kB][6]



### **[Metascape基本使用](https://metascape.org/gp/index.html#/main/step1)**

第一步：上传基因list，选择物种，并自定义功能注释

![image_1e6v925pag8vvo71lqphstd939.png-90.1kB][7]

第二步：选择功能富集选项

![image_1e6v95an71abd1lq51klf1omo29km.png-106.4kB][8]


第三步：选择富集的P value值，富集功能的选项等等，然后进行富集分析

![image_1e6vrh5irvtf1hvj1hfn5sigpj13.png-125.9kB][9]


第四步： 看结果，结果和DAVID的结果很类似

![image_1e6vs2vns1mld1ji3bni1lp81dk82a.png-193.7kB][10]



### **[ClusterProfiler使用](https://yulab-smu.github.io/clusterProfiler-book/)**

第一步： 下载和安装ClusterProfilre

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
```

第二步： 载入R包，进行单gene list的功能分析

```
library(clusterProfiler)
library(org.Hs.eg.db)
## 人的注释包，可以根据需要下载其它物种的注释包，User can query OrgDb online by AnnotationHub or build their own by AnnotationForge. An example can be found in the vignette of GOSemSim.

## ID 转换
RocA_up <- bitr(RocA_up$gene_name,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

## KEGG
kegg <- enrichKEGG(RocA_up$ENTREZID,organism = "hsa",keyType = "kegg",pvalueCutoff = 1,pAdjustMethod = "BH",qvalueCutoff = 0.5)
head(summary(kegg))
dotplot(kegg,font.size=10)

```
![image_1e6vspsf93v4f6i1m881i4ogrj2n.png-24.4kB][11]




第三步： 多基因list对比分析

```
## ID 转换
RocA_up <- bitr(RocA_up$gene_name,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
RocA_down <- bitr(RocA_down$gene_name,fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

## 构建比较矩阵
# up 
group_up <- rep("RUGs",times=length(RocA_up$ENTREZID))
data_up <- data.frame(geneID=RocA_up$ENTREZID,group=group_up)
# down
group_down <- rep("RDGs",times=length(RocA_down$ENTREZID))
data_down <- data.frame(geneID=RocA_down$ENTREZID,group=group_down)
# all
data_all_GO <- rbind(data_up,data_down)

## KEGG analysis
formula <- compareCluster(geneID~group, data=data_all_GO, fun='enrichKEGG',organism = "hsa",pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05)

## 绘图
y <- dotplot(formula,showCategory=20,includeAll=TRUE)
y

## 输出PPT
library(export)
graph2ppt(y,"KEGG_comparison.ppt")
```
![image_1e6vt9ekhhmm45p1m7q1oji1uoq34.png-28.9kB][12]



  [1]: http://static.zybuluo.com/sherking/na7yqaiy6q30pkc7nocsb50b/image_1e6v7iaqo1qg87lg14021etq1svn9.png
  [2]: http://static.zybuluo.com/sherking/9duqnmj6dtsht14t0oqs4s4o/image_1e6v7r2bc1ohuks1pdb1jn011l6m.png
  [3]: http://static.zybuluo.com/sherking/kklhhdxjpodl5uf3f3dic4qr/image_1e6v7tsbl1eeetk7kcd1vsv1b9e13.png
  [4]: http://static.zybuluo.com/sherking/uvgwrwgtmxnvwwf0k3ln0ln1/image_1e6v829pu35egbk1mqk13k7p2p1g.png
  [5]: http://static.zybuluo.com/sherking/vujg15yd9xjv07uk0m2owk37/image_1e6vdr6bt1ikv1urf1oet1kgs039.png
  [6]: http://static.zybuluo.com/sherking/8uku44vojellhjmll9kqvvy6/image_1e6vs1glf11ii15gocq41i422331t.png
  [7]: http://static.zybuluo.com/sherking/4dsldhhrd82niyv72snyb7tq/image_1e6v925pag8vvo71lqphstd939.png
  [8]: http://static.zybuluo.com/sherking/k7s5e59rdfrmwbd8wbcai4tx/image_1e6v95an71abd1lq51klf1omo29km.png
  [9]: http://static.zybuluo.com/sherking/2690he6uphc8qriuf9vy5ksb/image_1e6vrh5irvtf1hvj1hfn5sigpj13.png
  [10]: http://static.zybuluo.com/sherking/0llq7q641t8rnusi1ya534ow/image_1e6vs2vns1mld1ji3bni1lp81dk82a.png
  [11]: http://static.zybuluo.com/sherking/gyid341gd92bzb3ttfjn1v51/image_1e6vspsf93v4f6i1m881i4ogrj2n.png
  [12]: http://static.zybuluo.com/sherking/415qfjiluduvfu87mm9c1ay1/image_1e6vt9ekhhmm45p1m7q1oji1uoq34.png