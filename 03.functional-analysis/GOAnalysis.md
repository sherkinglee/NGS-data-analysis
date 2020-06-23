## **GO analysis**

Gene ontology (GO)分析在生物信息分析过程中非常常见，目前有很多方法可以进行GO分析，包括使用[DAVID](https://david.ncifcrf.gov/), [Metascape](https://metascape.org/gp/index.html#/main/step1)以及[ClusterProfiler](https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html)等等。这里我们分别进行介绍。

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

![image_1e6v8dqa31glggu163u368d5d2t.png-193kB][5]


第六步： 查看结果

![image_1e6v8jhu01pca1b9u1qhn8tpbve3t.png-95kB][6]



### **[Metascape基本使用](https://metascape.org/gp/index.html#/main/step1)**

第一步：上传基因list，选择物种，并自定义功能注释

![image_1e6v925pag8vvo71lqphstd939.png-90.1kB][7]

第二步：选择功能富集选项

![image_1e6v95an71abd1lq51klf1omo29km.png-106.4kB][8]


第三步：选择富集的P value值，富集功能的选项等等，然后进行富集分析

![image_1e6v9bj181t9bvc4ihb1go1qq29.png-152kB][9]


第四步： 看结果，结果和DAVID的结果很类似

![image_1e6v9qvmv1k4t2qf1bqm93tkj1m.png-212.2kB][10]



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

## BP分析
ego <- enrichGO(RocA_up$SYMBOL,OrgDb=org.Hs.eg.db,ont = "BP",keytype="SYMBOL",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.5)
dotplot(ego, font.size=10)

```
![image_1e6vc6g5jqoi3irakfpr55or9.png-32.1kB][11]

```
## CC分析
ego <- enrichGO(RocA_up$SYMBOL,OrgDb=org.Hs.eg.db,ont = "CC",keytype="SYMBOL",pvalueCutoff = 0.05,pAdjustMethod = "BH",qvalueCutoff = 0.5)
dotplot(ego, font.size=10)

```
![image_1e6vccv4ukg8vrpe7t1je7assm.png-27.3kB][12]


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

## GO analysis
formula <- compareCluster(geneID~group, data=data_all_GO, fun='enrichGO',OrgDb='org.Hs.eg.db',pvalueCutoff=0.05,pAdjustMethod = "BH",qvalueCutoff = 0.05,ont="BP")
summary(formula)
## 去除冗余terms
formula <- simplify(formula)

## 绘图
y <- dotplot(formula,showCategory=20,includeAll=TRUE)
y

## 输出PPT
library(export)
graph2ppt(y,"GO_comparison.ppt")
```
![image_1e6vd0ejv12tp1k616051cfg11b013.png-97.7kB][13]


### **ggplot2 富集分析**

这种情况适用于如果想利用Metascape或者DAVID的富集结果来画cluserProfiler的那种泡泡图，则可以考虑用这种脚本进行绘制。这个我还没有尝试过，以下脚本为旭东师兄所提供。

```
###########
###Type1###
###########
indata<-read.table("metascape_result_ggplot2.txt",header=T,sep="\t",check.names=F);
indata$Description <- factor(indata$Description, levels=as.character(indata$Description)); #reorder the factor level of Description
p = ggplot(indata,aes(y=Description,x=Rich_factor)) #ggplot2画图默认按照因子逆序从上到下画
pbubble = p + geom_point(aes(size=InList,color=-LogP))
pr = pbubble + scale_colour_gradient(low="yellow2",high="red") + 
  labs(color="-Log10(P-value)",size="Gene Number",x="Rich Factor",y="Terms Name",title="Top20 of genes enrichment functions") + 
   theme_bw()
pr

###########
###Type2###
###########
indata<-read.table("metascape_result_ggplot2.txt",header=T,sep="\t",check.names=F);
#indata$Description <- factor(indata$Description, levels=as.character(indata$Description));
###
p = ggplot(indata,aes(y=Description,x=InList))
pbubble = p + geom_point(aes(size=Rich_factor,color=-LogP))
pr = pbubble + scale_colour_gradient(low="yellow2",high="red") + 
  labs(color="-Log10(P-value)",size="Rich Factor",x="Gene Number",y="Terms Name",title="Top20 of genes enrichment functions") + 
   theme_bw()
pr

###########
###Type3###
###########
indata<-read.table("metascape_result_ggplot2.txt",header=T,sep="\t",check.names=F); #top20
#indata<-read.table("metascape_result_ggplot2.txt",header=T,sep="\t",check.names=F)[1:15,]; #top15
#indata<-read.table("metascape_result_ggplot2.txt",header=T,sep="\t",check.names=F)[1:10,]; #top10
###
p = ggplot(indata,aes(y=Description,x=0))
pbubble = p + geom_point(aes(size=Rich_factor+0.001,color=-LogP))
pr = pbubble + scale_colour_gradient(low="yellow2",high="red") + 
  labs(color="-Log10(P-value)",size="Rich Factor",x="Migrasome Enrichmed",y="Terms Name",title="Top20 of genes enrichment functions") + 
   theme_bw()
pr

```
  [1]: http://static.zybuluo.com/sherking/na7yqaiy6q30pkc7nocsb50b/image_1e6v7iaqo1qg87lg14021etq1svn9.png
  [2]: http://static.zybuluo.com/sherking/9duqnmj6dtsht14t0oqs4s4o/image_1e6v7r2bc1ohuks1pdb1jn011l6m.png
  [3]: http://static.zybuluo.com/sherking/kklhhdxjpodl5uf3f3dic4qr/image_1e6v7tsbl1eeetk7kcd1vsv1b9e13.png
  [4]: http://static.zybuluo.com/sherking/uvgwrwgtmxnvwwf0k3ln0ln1/image_1e6v829pu35egbk1mqk13k7p2p1g.png
  [5]: http://static.zybuluo.com/sherking/e9hbl13sjfm1idmneyomb808/image_1e6v8dqa31glggu163u368d5d2t.png
  [6]: http://static.zybuluo.com/sherking/qcovo4ugsgwrgx0inxvhw0zv/image_1e6v8jhu01pca1b9u1qhn8tpbve3t.png
  [7]: http://static.zybuluo.com/sherking/4dsldhhrd82niyv72snyb7tq/image_1e6v925pag8vvo71lqphstd939.png
  [8]: http://static.zybuluo.com/sherking/k7s5e59rdfrmwbd8wbcai4tx/image_1e6v95an71abd1lq51klf1omo29km.png
  [9]: http://static.zybuluo.com/sherking/v2ipqj98ddm7x4andvo3iofv/image_1e6v9bj181t9bvc4ihb1go1qq29.png
  [10]: http://static.zybuluo.com/sherking/6vi045jjtapw3bzpvsun8t2n/image_1e6v9qvmv1k4t2qf1bqm93tkj1m.png
  [11]: http://static.zybuluo.com/sherking/lwm9glr5zed5vrm9ieet2r4r/image_1e6vc6g5jqoi3irakfpr55or9.png
  [12]: http://static.zybuluo.com/sherking/v1478f8iqaixoyq4nq5hhqbo/image_1e6vccv4ukg8vrpe7t1je7assm.png
  [13]: http://static.zybuluo.com/sherking/utsipusg8x65c20c2mg9hc55/image_1e6vd0ejv12tp1k616051cfg11b013.png