## **venn plot**

韦恩图主要用来展示不同集合之间的重叠情况。我一般用R 语言来画韦恩图，主要用的package是VennDiagram，当然也可以选择其他的package，或者在线工具。

+ 安装VennDiagram

```
install.packages('VennDiagram')
```

+ 2D 韦恩图

```
# x1: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#x2: 1   2   3   4 200  34  67
#x3: 1  2  3  4  5  6  7  8  9 10 11 50 51 43 45 35 23 97 56 98
#x4: 1    2    3    4    5    6    7    8   51   43   45   35   23   97   56   98 1000  123 2353  545  634  232  254
#x5: 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  23  34  45  65  67 878  87  17  18  19

library(VennDiagram)

## 2D
venn.diagram(list(x1=x1,x2=x2),
             filename = "x1_x2.tiff",
             fill=c("red","blue"),
             alpha=c(0.5,0.5),
             cat.cex=1.5,
             cex=1.5,
             margin=c(0.1,0.1),
             cat.dist=c(0.1,0.1),
             cat.pos=c(-20,20),
             cat.fontfamily = "serif",
             cat.fontface = "bold")
```
![image_1ebi6mdad1u1f1vhvu7jrl8qf9.png-92.6kB][1]


+ 3D 韦恩图

```
# x1: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#x2: 1   2   3   4 200  34  67
#x3: 1  2  3  4  5  6  7  8  9 10 11 50 51 43 45 35 23 97 56 98
#x4: 1    2    3    4    5    6    7    8   51   43   45   35   23   97   56   98 1000  123 2353  545  634  232  254
#x5: 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  23  34  45  65  67 878  87  17  18  19

library(VennDiagram)

## 3D
venn.diagram(list(x1=x1,x2=x2,x3=x3),
             filename = "x1_x2_x3.tiff",
             fill=c("red","blue","green"),
             alpha=c(0.5,0.5,0.5),
             cat.cex=1.5,
             cex=1.5,
             margin=c(0.1,0.1),
             cat.col = c("darkred", "darkblue", "darkgreen"),
             cat.dist=c(0.1,0.1,0.1),
             cat.pos=c(-20,20,180),
             cat.fontfamily = "serif",
             cat.fontface = "bold")
```

![image_1ebi6p4gj2pstu21in01ksv1pkfm.png-130.2kB][2]

+ 3D 韦恩图

```
# x1: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#x2: 1   2   3   4 200  34  67
#x3: 1  2  3  4  5  6  7  8  9 10 11 50 51 43 45 35 23 97 56 98
#x4: 1    2    3    4    5    6    7    8   51   43   45   35   23   97   56   98 1000  123 2353  545  634  232  254
#x5: 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  23  34  45  65  67 878  87  17  18  19

library(VennDiagram)

## 4D

venn.diagram(list(x1=x1,x2=x2,x3=x3,x4=x4),
             filename = "x1_x2_x3_x4.tiff",
             col = "black",
             lty = "dotted", #边框线型改为"dotted"虚线
             lwd = 3, # 边框线的宽度
             fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
             alpha = 0.50,
             cex = 2.0,
             fontfamily = "serif",
             fontface = "bold",
             cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
             cat.cex = 1.8,
             cat.fontface = "bold",
             cat.fontfamily = "serif")
```

![image_1ebi6r732l2a1086shi1o7qra413.png-263.3kB][3]

+ 3D 韦恩图

```
# x1: 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#x2: 1   2   3   4 200  34  67
#x3: 1  2  3  4  5  6  7  8  9 10 11 50 51 43 45 35 23 97 56 98
#x4: 1    2    3    4    5    6    7    8   51   43   45   35   23   97   56   98 1000  123 2353  545  634  232  254
#x5: 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  23  34  45  65  67 878  87  17  18  19

library(VennDiagram)

venn.diagram(list(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5),
             filename = "x1_x2_x3_x4_x5.tiff",
             lty = "dotted",
             lwd = 2,
             col = "black",  #"transparent
             fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             alpha = 0.60,
             cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
             cat.cex = 1.5,
             margin = 0.07,
             cex = 1.5,
             cat.fontface = "bold",
             cat.fontfamily = "serif")
```

![image_1ebi6t8atcfj62ueh2jaga2u1g.png-265.7kB][4]

+ N D
超过5维数据集之间的比较，再用上述R包就不太适合了，这个时候UpSetR就显得十分强大。

```

## UpSetR
install.packages('UpSetR')
library(UpSetR)
library(RColorBrewer)
sets_bar_col=brewer.pal(length(test_data),"YlGn")
test_data <- fromList(list(x1,x2,x3,x4,x5)) ##转化为UpSetR需要的格式
colnames(test_data) <- paste0("x","",1:5) ##需要列名，否则会报错
x<-upset(test_data,nsets = 5,nintersects = 30,order.by = "freq",
      matrix.color = "#78C679", ## 矩阵中点的颜色
      main.bar.color =  topo.colors(12), ## bar plot的颜色
      mainbar.y.label="Intersection Size",
      sets.x.label="Gene Number",
      point.size=3, ##点大小
      line.size=1.5, ## 线粗细
      mb.ratio=c(0.6,0.4), ## 矩阵和barplot的比例
      text.scale=c(2.5,2.5,2.5,2.5,2.5,4),## ytitle, ylabel,xtitle,xlabel,sets,numbers
      sets.bar.color = sets_bar_col, ## size barplot 颜色
      shade.color="gray", ## 点的背景颜色
      set_size.show = T, ##是否展示gene number
      set_size.numbers_size=2) ## 数字大小
```

![image_1ebio6ct611ktg681anp17mfj399.png-67.1kB][5]


  [1]: http://static.zybuluo.com/sherking/d6iwageerudjxynrw9zu3lxt/image_1ebi6mdad1u1f1vhvu7jrl8qf9.png
  [2]: http://static.zybuluo.com/sherking/mamdolqm4nykzwcuushmnuyt/image_1ebi6p4gj2pstu21in01ksv1pkfm.png
  [3]: http://static.zybuluo.com/sherking/1gw7vkroo8fyvdyoom33ck5e/image_1ebi6r732l2a1086shi1o7qra413.png
  [4]: http://static.zybuluo.com/sherking/c36n52ee3tmhozrakcpeilt6/image_1ebi6t8atcfj62ueh2jaga2u1g.png
  [5]: http://static.zybuluo.com/sherking/v6lrzwm8yym7gokppinesx87/image_1ebio6ct611ktg681anp17mfj399.png