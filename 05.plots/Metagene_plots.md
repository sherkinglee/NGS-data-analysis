## **Metagene plot**

Metagene plot是指把所有转录本或者基因按照某种规则（如按照某个特定的位置）对齐，画出每个位置上的相对度量。比如说，在ribosome profiling数据中想要看每个转录本起始密码子附近的RPF的相对表达值，就可以用metagene plot的形式展现出来。这里可以提供R和Python两种绘图方式。

### **Python 绘制metagene plot**

关于这一块，主要还是针对Ribo-seq数据来说的，因为我博士阶段主要分析的数据类型就是Ribo-seq，计算的工具主要也是自己开发的小程序[RiboMiner](https://github.com/xryanglab/RiboMiner).

+ 数据准备阶段,具体产生方式见[RiboMiner](https://github.com/xryanglab/RiboMiner)或者[CodeOcean::Ribomier](https://codeocean.com/capsule/3780896/tree/v1)

```
## 数据主要是RiboMiner生成的，如RiboMiner::MetageneAnalysis输出的结果
sample  start_density   stop_density
si-Ctrl-1       1.0235216043853657      0.21732726629052093
si-Ctrl-1       0.6689095490141639      0.3019296172715671 
...
...
si-eIF5A-2      0.1802313319028725      0.6756468494549336 
si-eIF5A-2      0.1659164008524256      0.6206370659429937
```

+ python数据绘图，具体脚本内容见[RiboMimer::PlotMetageneAnalysis](https://github.com/xryanglab/RiboMiner/blob/master/RiboMiner/PlotMetageneAnalysis.py)

```
PlotMetageneAnalysis -i $results/MetageneAnalysis/eIF5A_CDS_normed_dataframe.txt -o $results/MetageneAnalysis/eIF5A_CDS_normed -g si-Ctrl,si-eIF5A -r si-Ctrl-1,si-Ctrl-2__si-eIF5A-1,si-eIF5A-2 -U nt -u 50 -d 50 --mode mean --CI 0.95 --axhline 1
```

+ python绘图结果展示

![image_1ebfjl6mep7m1suv9g6fdc14iu9.png-225.7kB][1]



### **R 绘制metagene plot**
输入的数据都是RiboMiner生成的文件，只不过方便绘图，才用R语言进行绘制了一次。

```
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))

## 读数据
all_density <- read.table("paper_RPF_mean_density_dataframe.txt",sep="\t",header = T,stringsAsFactors = F)
Stop_codon_density <- read.table("paper_RPF_stop_density_dataframe.txt",header = T,sep="\t",stringsAsFactors = F)
UTR_density <- read.table("paper_RPF_utr_density_dataframe.txt",header=T,sep="\t",stringsAsFactors = F)
motif_density <- read.table("paper_RPF_motif_motifDensity_dataframe.txt",sep="\t",header = T,stringsAsFactors = F)

dim(all_density)
head(all_density)

id<-rep(seq(0,500,1),times=2)
all_density <-mutate(all_density,id=id)

g <- ggplot(data=all_density,aes(x=id,y=start_density,color=sample)) +geom_line() + labs(x="Distance from start codon",y="Relative footprint density (AU)",title="Ribosome footprint density profiles") + theme_bw()+theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks = seq(0,500,50),labels = seq(0,500,50))
g

```

结果展示如下：

![image_1ebfk3mnf1qk917fu15t91kk2s03m.png-122.8kB][2]


两组数据比较的话，可以用下面的方式实现

```
control_1 <- Stop_codon_density[Stop_codon_density$sample=="SRR5008134_wt",]
control_2 <- Stop_codon_density[Stop_codon_density$sample=="SRR5008135_wt",]
treat_1 <- Stop_codon_density[Stop_codon_density$sample=="SRR5008136_d",]
treat_2 <- Stop_codon_density[Stop_codon_density$sample=="SRR5008137_d",]

control_mean <- data.frame(sample=rep("control",nrow(control_1)),start_density=2.5*(control_1$start_density+control_2$start_density)/2,stop_density=2.5*(control_1$stop_density+control_2$stop_density)/2)

treat_mean <- data.frame(sample=rep("eIF5Ad",nrow(treat_1)),start_density=2.5*(treat_1$start_density+treat_2$start_density)/2,stop_density=2.5*(treat_1$stop_density+treat_2$stop_density)/2)

Stop_codon_density_mean <- rbind(control_mean,treat_mean)
id <- rep(seq(-100,50,length.out = nrow(control_mean)),times=2)

Stop_codon_density_mean <- mutate(Stop_codon_density_mean,id=id)

p <-  ggplot(data=Stop_codon_density_mean,aes(x=id,y=stop_density,color=sample)) +geom_line() + facet_grid(sample~.)+labs(x="Distance from stop codon (nt)",y="Relative footprint density (AU)",title="Ribosome footprint density profiles") + theme_bw()+theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks = seq(-100,50,50),labels = seq(-100,50,50))
p 
```
![image_1ebfk8fjn11gr8c7tn7ijr1gsg13.png-164.2kB][3]


对于特定motifs附近的density情况也可以用下面方式实现：

```
motif_density <- mutate(motif_density,control_mean=(SRR5008134_wt+SRR5008135_wt)/2,eIF5Ad_mean=(SRR5008136_d+SRR5008137_d)/2)
PP <- motif_density[motif_density$motif=="PP",]
PPP <- motif_density[motif_density$motif=="PPP",]

PP_shaped <- melt(PP,id.vars =colnames(PP)[1:5])
PPP_shaped <- melt(PPP,id.vars = colnames(PPP)[1:5])
id <- rep(seq(-50,50,length.out = 101),times=2)

PP_shaped <- mutate(PP_shaped,id=id)
PPP_shaped <- mutate(PPP_shaped,id=id)
p1 <- ggplot(data=PP_shaped,aes(x=id,y=value,color=variable)) + geom_line()
p1 + ylim(0,0.15) + labs(x="Distance of P-site to PP (nt)",y="Relative footprint density",title="Ribosome footprint density") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +guides(color=guide_legend(title=NULL))

p2 <- ggplot(data=PPP_shaped,aes(x=id,y=value,color=variable)) + geom_line()
p2 + ylim(0,0.15) + labs(x="Distance of P-site to PPP (nt)",y="Relative footprint density",title="Ribosome footprint density") + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +guides(color=guide_legend(title=NULL))
```

![image_1ebfkfntt36h1hpr1qe1bij1ogp1g.png-50.8kB][4]


  [1]: http://static.zybuluo.com/sherking/0bi66vxjf97wax8itxiu055v/image_1ebfjl6mep7m1suv9g6fdc14iu9.png
  [2]: http://static.zybuluo.com/sherking/h2nz4094jslkxuzxs0ivfff3/image_1ebfk3mnf1qk917fu15t91kk2s03m.png
  [3]: http://static.zybuluo.com/sherking/c2gvb6rp88l7pjrsoqlv1u3p/image_1ebfk8fjn11gr8c7tn7ijr1gsg13.png
  [4]: http://static.zybuluo.com/sherking/b1wjqpbyzy52515dnwr2eaoi/image_1ebfkfntt36h1hpr1qe1bij1ogp1g.png