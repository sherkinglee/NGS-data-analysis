## **将Rplots导出为PPT格式**

我们平时用R语言，python或者其他工具绘制的图表一般都是用作课题讨论展示时所用，真正到了文章发表的时候，所有程序绘制的图表为了符合发表要求都需要重新进行绘制。然而，用脚本调参绘制图表非常繁琐，因此能够对图表每个细节进行修改的office图表就显得很方便。当然有人会argue说用AI调图也很方便啊，但是AI这个付费软件不是每个人都用得起的，即使有破解版，也耗很大的内存，而且也不是每个人都精通AI图表绘制。相对而言，office技巧就显得大众许多。我一开始也不知道R语言的图还能转化为PPT，直到实验室微沙龙的时候师姐分享了一个R包export可以将Rplots导出为pptx格式，随后我尝试了一下感觉还不错，但是export似乎只能导出ggplot2等画出的图表，而基础包画的图导不出来。而且因为作者没有持续更新，export目前也被CRAN剔除掉了。后面又看到Y叔的公众号，分享了eoffice包也可以实现类似的功能。所以这里分享这两个包的使用。

### **export使用**

```
install.packages('export') ## 目前下载不了，因为被CRAN剔除了（执笔时20200624）
library(export)
graph2ppt(p,"volcano_plot.ppt") # p为ggplot等产生的对象

```


### **eoffice使用**

```
## export被下架，所以可以选择另外一个eoffice,结合ggplotify可以把任何R plot输出为PPT格式
install.package('eoffice')
library(eoffice)
topptx(g,'volvano_plot2.pptx')

## 如果是基础包画出的图形，可以用ggplotify转化为ggplot类型图标
g=as.ggplot(plot(x,y))
topptx(g,'scatter_plot.pptx')

```