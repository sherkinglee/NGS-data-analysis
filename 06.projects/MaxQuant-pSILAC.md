# MaxQuant 分析pSILAC数据

标签：蛋白质组学 bioinformatics pSILAC MaxQuant 20211006

[TOC]

## 概述

蛋白质组学数据可以用来定量不同条件下的蛋白质变化情况。目前有很多分析蛋白组学数据的工具或者软件（详见[蛋白组学数据分析工具](https://www.cnblogs.com/jessepeng/p/13575757.html))，包括[MaxQuant](https://www.maxquant.org/), [PD](https://www.cnblogs.com/jessepeng/p/13579947.html)等等。在分析小鼠早期胚胎发育过程中的翻译组学数据时，涉及到部分蛋白质组学定量的分析。这里总结下MaxQuant分析pSILAC质谱数据的方法和过程。

## 数据来源

这里所用的数据比较多，是2019年发表在scientific report上的一篇工作。文章链接：[Israel, S., et al. (2019). An integrated genome-wide multi-omics analysis of gene expression dynamics in the preimplantation mouse embryo. Sci Rep-Uk 9.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6746714/).数据存储在[ProteomeXchange](http://www.proteomexchange.org/), accession number 为[PXD007082](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD007082)。获取[数据表格](https://github.com/sherkinglee/NGS-data-analysis/tree/master/06.projects/Download.txt)，准备下载。

```
## Download.txt
ID      NAME    URI     TYPE    MAPPINGS
9       20140227_exp0616_Michele_embryo-F9_HPLC_0_01.raw
10      20140227_exp0616_Michele_embryo-F9_HPLC_0_02.raw
11      20140227_exp0616_Michele_embryo-F9_HPLC_0_03.raw
12      20140227_exp0616_Michele_embryo-F9_HPLC_0_04.raw
13      20140227_exp0616_Michele_embryo-F9_HPLC_0_05.raw
14      20140227_exp0616_Michele_embryo-F9_HPLC_0_06.raw
15      20140227_exp0616_Michele_embryo-F9_HPLC_0_07.raw
16      20140227_exp0616_Michele_embryo-F9_HPLC_0_08.raw
17      20140227_exp0616_Michele_embryo-F9_HPLC_0_09.raw
18      20140227_exp0616_Michele_embryo-F9_HPLC_0_10.raw
19      20140227_exp0616_Michele_embryo-F9_HPLC_0_11.raw
20      20140227_exp0616_Michele_embryo-F9_HPLC_0_12.raw
21      20140227_exp0616_Michele_embryo-F9_HPLC_0_13.raw
22      20140227_exp0616_Michele_embryo-F9_HPLC_0_14.raw
23      20140227_exp0616_Michele_embryo-F9_HPLC_0_15.raw

less -S Download.txt|while read ID NAME URI TYPE MAPPINGs;do

URL="http://"${URI:6}
#echo $URL
nohup wget -c $URL &
done

```

## [MaxQuant安装](https://www.maxquant.org/)

+ [MaxQuant本地安装](https://www.maxquant.org/)需要注册，然后通过邮件发送的链接，进行下载安装。
+ [MaxQuant服务器安装](https://anaconda.org/bioconda/maxquant/files)可以通过conda进行安装。

```
conda install -c bioconda maxquant
```
这里我们直接通过下载MaxQuant，然后传到服务器进行使用

```
wget -c https://anaconda.org/bioconda/maxquant/2.0.1.0/download/noarch/maxquant-2.0.1.0-py39hdfd78af_2.tar.bz2

tar -xvjf maxquant-2.0.1.0-py39hdfd78af_2.tar.bz2
chmod 755 bin/*
```
+ 安装完MaxQuant之后，在Linux使用，还需要下载[Mono](https://www.mono-project.com/)，这里我们还是通过conda安装

```
conda install mono

mono -h
```

+ 安装完mono之后，在终端测试MaxQuant是否能够正常使用：

```
mono ./bin/MaxQuantCmd.exe -h

MaxQuantCmd 2.0.1.0
Copyright © Max-Planck-Institute of Biochemistry 2021

ERROR(S):
  Option 'h' is unknown.
  A required value not bound to option name is missing.
USAGE:
Complete run of an existing mqpar.xml file:
  MaxQuantCmd.exe mqpar.xml
Print job ids/names table:
  MaxQuantCmd.exe --dryrun mqpar.xml
Rerunning second and third parts of the analysis:
  MaxQuantCmd.exe --partial-processing-end=3
...
```

+ 如果出现.NET缺失，提示下载.NET的时候，还是通过conda下载安装核心部件

```
conda install -c conda-forge dotnet-sdk
```

有时候上述方法还可能下载不下来，所以直接通过官网下载[dotnet x64](https://dotnet.microsoft.com/download/dotnet/2.1).

```
tar  -xvzf dotnet-sdk-2.1.818-linux-x64.tar.gz
./dotnet -h
```

当上述过程均没问题，mono, dotnet, MaxQuant全部下载之后，开始准备运行MaxQuant.

## MaxQuant本地使用

+ 通过邮件链接下载MaxQuant之后，直接解压缩，运行.exe文件就可以打开MaxQuant

![image_1fhak0jn41r1ltqrr0m15ia1rc79.png-37.2kB][1]

+ 点击“原始数据”下的load即可以上传原始.raw格式文件。

![image_1fhak3vms10ao1jfh11p31mthdt9m.png-147.1kB][2]

+ 如果数据比较少，通过set experiment以及set fraction等来设置实验和fraction。parameter group基本不用变，同一个样本设置同一个experiment名字，一个experiment有不同的fraction分别设置对应的fraction编号，如这套数据每个样本是20个fraction，fraction设置为1-20.
+ 如果数据比较多，如这套数据一个replicate有7个时期的样本，包括M2, 1cell, 2cell, 4cell, 8cell, 16cell(morula),32cell(blastocyst)，每个时期20个fraction,手动设置比较复杂，可以提前在txt文件中设置好相关内容，如：
![image_1fhakdr5u1ud51hveg7iintu4k1g.png-190.8kB][3]

+ 然后点击read from file,即可：
![image_1fhakg1kgqbg1tbgucf167n1b011t.png-156.7kB][4]

+ 然后点击“组特定参数”，Type设置情况如下：standard, 2,heavy label 设置Arg10, lys8
![image_1fhakiv3ud0r1kvq1hhef4326g2a.png-25.5kB][5]

modification参数默认
![image_1fhakl3gqp94125gff1gi910j32n.png-34.5kB][6]

![image_1fhakltvc193tlq8nhg1tog1q9m34.png-11.7kB][7]
![image_1fhakmsb3hp4um91ivfon61oas3h.png-25.4kB][8]

+ 点击“全局参数”，sequence添加数据库文件，从uniprot数据库下载

![image_1fhakp756d1u1df6od0appgqp3u.png-34.4kB][9]

然后identification下勾选match between run

![image_1fhakrcto1kdu10fhetg1qav15ta4b.png-27.6kB][10]

+ 然后点击“性能”show activities

![image_1fhakss121tv8fpup6othl1anh4o.png-13.8kB][11]

+ 最后点击开始，既可以运行MaxQuant

![image_1fhakva4ktv31l3kqa0ej11e5355.png-220.1kB][12]

+ 在点击开始之后，就可以产生一个mqpar.xml文件，这里存储maxquant运行的参数，数据库文件，raw文件等信息，可以将之转移到Linux服务器，替换数据库地址，raw文件地址之后在linux下运行MaxQuant。


## MaxQuant 服务器使用

+ 安装之前的方式，安装好MaxQuant，这里我们使用的服务器MaxQuant版本是maxquant-2.0.1.0，dotnet 版本是2.1.818
+ 将本地运行MaxQuant产生的mqpar.xml文件上传服务器，然后更改文件内部数据库和raw文件地址，在vim中用%s全局替换即可。
+ 运行MaxQuant

```
mono /workdata/home/lifj/lifj/software/MaxQuant/maxquant-2.0.1.0/bin/MaxQuantCmd.exe  mqpar.xml
```

+ 查看结果,显示运行流程。速度比本地块，所以大数据不建议在本地跑MaxQuant，在本地用MaxQuant产生mqpar.xml文件之后，传到linux上运行即可。
```
$ ls combined/proc/|less -S
Assemble_run_info 001140.finished.txt
Assemble_run_info 002140.finished.txt
Assemble_run_info 003140.finished.txt
Assemble_run_info 004140.finished.txt
Assemble_run_info 005140.finished.txt
Assemble_run_info 006140.finished.txt
Assemble_run_info 007140.finished.txt
Assemble_run_info 008140.finished.txt
Assemble_run_info 009140.finished.txt
Assemble_run_info 010140.finished.txt
Assemble_run_info 011140.finished.txt
Assemble_run_info 012140.finished.txt
Assemble_run_info 013140.finished.txt
Assemble_run_info 014140.finished.txt
Assemble_run_info 015140.finished.txt
Assemble_run_info 016140.finished.txt

```

  [1]: http://static.zybuluo.com/sherking/ymqm8cur8et4ay3zwro700g5/image_1fhak0jn41r1ltqrr0m15ia1rc79.png
  [2]: http://static.zybuluo.com/sherking/abtwpqabqtkeusytxnsrhfnu/image_1fhak3vms10ao1jfh11p31mthdt9m.png
  [3]: http://static.zybuluo.com/sherking/icggtjyh658ydilckvarrgc5/image_1fhakdr5u1ud51hveg7iintu4k1g.png
  [4]: http://static.zybuluo.com/sherking/mlbrrmuz3qsbtv1vn9i75m4d/image_1fhakg1kgqbg1tbgucf167n1b011t.png
  [5]: http://static.zybuluo.com/sherking/nh3goqr31xjdv4slk7ipqcf9/image_1fhakiv3ud0r1kvq1hhef4326g2a.png
  [6]: http://static.zybuluo.com/sherking/bq85z0kdasxnlj1nfnju8wwf/image_1fhakl3gqp94125gff1gi910j32n.png
  [7]: http://static.zybuluo.com/sherking/z5boh5hnqwuptxg2ekw7g1j9/image_1fhakltvc193tlq8nhg1tog1q9m34.png
  [8]: http://static.zybuluo.com/sherking/7wgqbfjl2zo9qv1lzxo8b5s6/image_1fhakmsb3hp4um91ivfon61oas3h.png
  [9]: http://static.zybuluo.com/sherking/jan6gdpw1t7r04i2sipmh72d/image_1fhakp756d1u1df6od0appgqp3u.png
  [10]: http://static.zybuluo.com/sherking/637uj3zdoxaltbaqqfzp0xqv/image_1fhakrcto1kdu10fhetg1qav15ta4b.png
  [11]: http://static.zybuluo.com/sherking/46h10bi42hu218nbll07ggwc/image_1fhakss121tv8fpup6othl1anh4o.png
  [12]: http://static.zybuluo.com/sherking/vttl9tcfpfky9or3oh2edwby/image_1fhakva4ktv31l3kqa0ej11e5355.png