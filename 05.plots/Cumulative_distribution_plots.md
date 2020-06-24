## **累积密度分布图： cumulative distribution plot**

累积密度分布图可以有多种方法进行绘制，这里我们介绍三种： R+sROF, R+ggpubr,python。

+ 数据准备

```
## 数据类型，dataframe，大致长这样
           samples FC(RocA003/DMSO) FC(RocA03/DMSO) FC(RocA3/DMSO)
1 RocA-insensitive        0.8129606        0.000000       3.293892
3 RocA-insensitive        5.8550586       10.831858      20.100356
5 RocA-insensitive        1.2461415        1.864746       1.238724
6 RocA-insensitive        0.0000000        0.000000       1.176471
8 RocA-insensitive        3.5696203        3.686275       1.035813
```

+ R+sROC

```
library("sROC");
plot(kCDF(RPFdist[,2])$x,kCDF(RPFdist[,2])$Fhat,type="l",lwd=3,col="red",xlab="FC(RocA/DMSO)",ylab="Cumulative Fraction",main="Rep1",frame.plot=FALSE)
lines(kCDF(RPFdist[,3])$x,kCDF(RPFdist[,3])$Fhat,lwd=3,col="black")
lines(kCDF(RPFdist[,4])$x,kCDF(RPFdist[,4])$Fhat,lwd=3,col="green")
legend(40,0.4,lty=c(1,1),legend=c("RocA003","RocA03","RocA3"),col=c("red","black","green"),bg="aliceblue",lwd=2)
```

![image_1ebiu5qes181f2l22de1ftg1hln9.png-47.5kB][1]


+ R+ggpubr

```
## 数据准备
RPFdist_FC <- RPFdist[,c(1,6,7,8)]
RPFdist_raw <- RPFdist[,c(1,2,3,4,5)]
RPFdist_melt <- melt(RPFdist_FC,id.vars = "samples")

## RPFdist_raw
           samples        DMSO    RocA003     RocA03      RocA3
1 RocA-insensitive 0.026171741 0.02127660 0.00000000 0.08620690
3 RocA-insensitive 0.009232026 0.05405405 0.10000000 0.18556701
5 RocA-insensitive 0.043247269 0.05389222 0.08064516 0.05357143
6 RocA-insensitive 0.050000000 0.00000000 0.00000000 0.05882353
8 RocA-insensitive 0.015957447 0.05696202 0.05882353 0.01652893
9 RocA-insensitive 0.012531328 0.00000000 0.00000000 0.08620690

## RPFdist_FC
           samples FC(RocA003/DMSO) FC(RocA03/DMSO) FC(RocA3/DMSO)
1 RocA-insensitive        0.8129606        0.000000       3.293892
3 RocA-insensitive        5.8550586       10.831858      20.100356
5 RocA-insensitive        1.2461415        1.864746       1.238724
6 RocA-insensitive        0.0000000        0.000000       1.176471
8 RocA-insensitive        3.5696203        3.686275       1.035813
9 RocA-insensitive        0.0000000        0.000000       6.879310
## RPFdist_melt
           samples         variable     value
1 RocA-insensitive FC(RocA003/DMSO) 0.8129606
2 RocA-insensitive FC(RocA003/DMSO) 5.8550586
3 RocA-insensitive FC(RocA003/DMSO) 1.2461415
4 RocA-insensitive FC(RocA003/DMSO) 0.0000000
5 RocA-insensitive FC(RocA003/DMSO) 3.5696203
6 RocA-insensitive FC(RocA003/DMSO) 0.0000000

```

```
## 绘制图片
RPFdist_melt <- melt(RPFdist_raw,id.vars = "samples")
p<-ggecdf(RPFdist_melt,x="value",color="samples",palette = "npg",font.label = 16,size = 2) + labs(x="log2(RPFdist)",y="Cumulative fraction") 
pp <- ggpar(p,legend = "right",legend.title = "Samples",xscale="log2")
pp
```
![image_1ebiuhmtr1n1g1o8i95t1d92alnm.png-98.1kB][2]

```
## 输出到PPT
install.package('eoffice')
library(eoffice)
topptx(pp,'ecdf.pptx')
```

+ python

```
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from scipy import stats

## prepare data for cdf plot
hist1,bins1=np.histogram(up_data.iloc[:,-1],1000,density=True)
cdf1=np.cumsum(hist1*np.diff(bins1))
hist2,bins2=np.histogram(up_mt_data.iloc[:,-1],1000,density=True)
cdf2=np.cumsum(hist2*np.diff(bins2))
hist3,bins3=np.histogram(down_data.iloc[:,-1],1000,density=True)
cdf3=np.cumsum(hist3*np.diff(bins3))
hist4,bins4=np.histogram(unblocked_data.iloc[:,-1],1000,density=True)
cdf4=np.cumsum(hist4*np.diff(bins4))

coor1=stats.mannwhitneyu(up_data.iloc[:,-1],unblocked_data.iloc[:,-1]) # un-up
coor2=stats.mannwhitneyu(up_mt_data.iloc[:,-1],unblocked_data.iloc[:,-1]) # un-up_mt
coor3=stats.mannwhitneyu(down_data.iloc[:,-1],unblocked_data.iloc[:,-1]) # un-down
coor4=stats.mannwhitneyu(down_data.iloc[:,-1],up_data.iloc[:,-1]) # up-down
coor5=stats.mannwhitneyu(up_mt_data.iloc[:,-1],up_data.iloc[:,-1]) # up-up_mt

## plot cdf
plt.rc("font",weight="bold")
fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(bins1[:-1],cdf1,color='r',label='up',linewidth=1.5)
# ax.plot(bins2[:-1],cdf2,color='b',label='up_mt')
ax.plot(bins3[:-1],cdf3,color='g',label="down",linewidth=1.5)
ax.plot(bins4[:-1],cdf4,color='k',label='unblocked',linewidth=1.5)
ax.text(0.8,0.4,r"P(up_and_unblocked)="+str('%2g'%coor1[1]),fontdict={'weight':'bold',"size":12,"family":"Arial"})
ax.text(0.8,0.3,r"P(down_and_unblocked)="+str('%2g'%coor3[1]),fontdict={'weight':'bold',"size":12,"family":"Arial"})
ax.text(0.8,0.2,r"P(up_and_down)="+str('%2g'%coor4[1]),fontdict={'weight':'bold',"size":12,"family":"Arial"})
ax.set_xlabel("log2FC(si-eIF3ed/si-control)",fontdict={'weight':'bold',"size":20,"family":"Arial"})
ax.set_ylabel("Cummulative Fraction",fontdict={'weight':'bold',"size":20,"family":"Arial"})
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.tick_params(which="both",width=2,labelsize=12)
plt.legend(loc="best",prop={'weight':'bold',"size":20,"family":"Arial"})
plt.show()

```

![image_1ebivii4e9a11tl910qac0m1kcj9.png-75.4kB][3]


  [1]: http://static.zybuluo.com/sherking/812pgo27isf5jhzxhbo4s3wh/image_1ebiu5qes181f2l22de1ftg1hln9.png
  [2]: http://static.zybuluo.com/sherking/1p2gm9va0c1ru4a3lvn01jhz/image_1ebiuhmtr1n1g1o8i95t1d92alnm.png
  [3]: http://static.zybuluo.com/sherking/lgss0koyfv6gto3qtsoprcjf/image_1ebivii4e9a11tl910qac0m1kcj9.png