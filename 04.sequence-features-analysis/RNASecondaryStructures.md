## **RNA secondary structures**

这部分我也没有太多的经验，主要是用[RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)计算过5'UTR的最小自由能，统计过不同转录本形成二级结构可能性的metagene结果等等。

+ [RNAfold Download](https://www.tbi.univie.ac.at/RNA/index.html) and installation
+ 单个sequence的二级结构预测

```
RNAfold  -i sequence.fa -p --outfile
```

+ 多个不同序列的二级结构预测

```
#!/usr/bin/python
# -*- coding:UTF-8 -*-

import os
import sys
import pandas as pd
import numpy as np
from collections import defaultdict
from itertools import groupby
from optparse import OptionParser



def creat_parse():
    usage="python %prog [options]"
    parser=OptionParser(usage=usage)
    parser.add_option("-i","--input",action="store",type='string',dest="input",help="Input file is a directoty where the input files are s
    parser.add_option("-o","--output",action="store",type="string",dest="output",help="Output is a directory where the results will be sto
    return parser

def CallRNAfold(input,output):
    if not os.path.exists(output):
        os.makedirs(output)
    files=os.listdir(input)
    files=[file for file in files if file.endswith(".fasta")]
    os.chdir(output)
    for file in files:
        cmd="RNAfold" + ' -i '+'../fastaFiles'.strip()+"/"+file + ' -p --outfile '
        os.system(cmd)
def main():
    parser=creat_parse()
    (options,args)=parser.parse_args()
    inputDir,outputDir=options.input,options.output
    print("Start rnafold...",file=sys.stderr)
    CallRNAfold(inputDir,outputDir)
    print("Finish rnafold...",file=sys.stderr)

if __name__=="__main__":
    main()

## 脚本用法
python CallRNAfold.py -i inputDir -o outputDir

## 其中inputDir包括每条序列的fasta文件，outputDir包括每条序列输出的二级结构结果
```


+ 提取每个序列的最小自由能（MFE）

```
#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd
import numpy as np
import os
import sys
import re
from optparse import OptionParser


def creat_parse():
    usage="python %prog [options]"
    parser=OptionParser(usage=usage)
    parser.add_option("-i","--input",action="store",type='string',dest="input",help="Input file is a directoty where the input files are s
tored.")
    parser.add_option("-o","--output",action="store",type="string",dest="outPrefix",help="Output prefix of the final results")
    return parser

def preprocess(inputDir,outPrefix):
    fileList=os.listdir(inputDir)
    files=[x for x in fileList if x.endswith(".fold")]
    fout=open(outPrefix+"_MEF.txt",'w')
    fout.write("%s\t%s\n" %("transcript","MEF"))
    for fn in files:
        geneName=fn.split(".fold")[0]
        #print(geneName)
        with open(inputDir.rstrip("/")+"/"+fn,'r') as f:
            lines=f.readlines()
            i=0
            for line in lines:
                i+=1
                if i==3:
                    mef=re.findall("(-?\d+.\d+)",line)[0]
                    fout.write("%s\t%s\n" % (geneName,str(mef)))
    fout.close()

def main():
    parser=creat_parse()
    (options,args)=parser.parse_args()
    (input,outPrefix)=(options.input,options.outPrefix)
    print("Start preprocessing the data...",file=sys.stderr)
    preprocess(input,outPrefix)
    print("Finish the step of processing!")

if __name__=="__main__":
    main()

## 脚本使用
python extractMFEFromRNAfold.py -i inputDir -o outprefix

##其中inputDir则为上述脚本生成的，-o为输出文件的前缀
## 输出结果为
transcript      MEF
ENST00000270460 -163.10
ENST00000537403 -306.50
ENST00000361752 -331.50
ENST00000256079 -150.80
ENST00000367106 -188.00
ENST00000374403 -21.80
ENST00000361436 -32.70

```

+ 提取每条序列在每个位置上形成二级结构的可能性

```
#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd
import numpy as np
import os
import sys
from optparse import OptionParser


def creat_parse():
    usage="python %prog [options]"
    parser=OptionParser(usage=usage)
    parser.add_option("-i","--input",action="store",type='string',dest="input",help="Input file is a directoty where the input files are stored.")
    parser.add_option("-o","--output",action="store",type="string",dest="outPrefix",help="Output prefix of the final results")
    parser.add_option("-f","--fastaFile",action="store",type="string",dest="fastaFile",help="fasta file containing all the transcript sequences we wanna analyze")
    return parser


def getLengthDict(fastaFile):
    transcriptDict={}
    with open(fastaFile,'r') as f:
        lines=f.readlines()
        for i in np.arange(len(lines)):
            if ">" in lines[i]:
                geneName=lines[i].strip(">").strip().split(" ")[0]
                sequence=lines[i+1].strip()
                length=len(sequence)
                transcriptDict[geneName]=length
    return transcriptDict

def preprocess(inputDir,outPrefix,fastaFile):
    fileList=os.listdir(inputDir)
    files=[x for x in fileList if x.endswith("dp.ps")]
    templbox,tempubox=np.zeros([1,3]),np.zeros([1,3])
    transcriptDict=getLengthDict(fastaFile)
    for fn in files:
        geneName=fn.split("_")[0]
        gene_length=transcriptDict[geneName]
        if gene_length < 0: ## adjust based on your need
            continue
        #print("gene length: ",gene_length)
        lbox,ubox=[],[]
        with open(inputDir.rstrip("/")+"/"+fn,'r') as f:
            lines=f.readlines()
            for line in lines:
                tmp=line.strip().split(" ")
                if len(tmp)==4 and tmp[3]=="ubox":
                    pos,target,prob=int(tmp[0]),int(tmp[1]),float(tmp[2])**2
                    ubox.append([pos,target,prob])
                if len(tmp)==4 and tmp[3]=="lbox":
                    pos,target,prob=int(tmp[0]),int(tmp[1]),float(tmp[2])
                    lbox.append([pos,target,prob])
        lbox=pd.DataFrame(lbox,columns=["pos","target","probability"])
        #print("lobx: ",lbox.shape)
        lbox.head(10)
        ubox=pd.DataFrame(ubox,columns=["pos","target","probability"])
        #print("uobx: ",ubox.shape)
        ubox.head(10)
        lbox_grouped=lbox.groupby("pos").max().reset_index()
        #print("lbox_grouped:",lbox_grouped.shape)
        lbox_grouped.head(10)
        ubox_grouped=ubox.groupby("pos").max().reset_index()
        #print("rbox_grouped:",ubox_grouped.shape)

        lbox=lbox_grouped[['pos','probability']].fillna(value=0)
        lbox.head(10)
        ubox=ubox_grouped[['pos','probability']].fillna(value=0)
        ubox.head(10)

        pos_list=[[i] for i in np.arange(gene_length)]
        pos_data=pd.DataFrame(pos_list,columns=['pos'])
        #print("pos_data: ",pos_data.shape)
        pos_data.head(10)
        ubox_data=pd.merge(pos_data,ubox,how="left").fillna(value=0)
        lbox_data=pd.merge(pos_data,lbox,how='left').fillna(value=0)
        ubox_data=ubox_data.sort_values("pos")
        lbox_data=lbox_data.sort_values("pos")


        geneList_l=((geneName+" ")*lbox_data.shape[0]).strip().split(" ")
        geneList_u=((geneName+" ")*ubox_data.shape[0]).strip().split(" ")
        geneList_l=pd.DataFrame(np.array(geneList_l).reshape(-1,1),columns=['gene'])
        geneList_u=pd.DataFrame(np.array(geneList_u).reshape(-1,1),columns=['gene'])

        lbox_data=pd.merge(geneList_l,lbox_data,left_index=True,right_index=True)
        ubox_data=pd.merge(geneList_u,ubox_data,left_index=True,right_index=True)
        #print("lbox_data: ",lbox_data.shape)
        #print("ubox_data: ",ubox_data.shape)

        templbox=np.vstack((templbox,lbox_data))
        tempubox=np.vstack((tempubox,ubox_data))

    lbox_probs=pd.DataFrame(templbox,columns=["gene","pos","probability"])
    ubox_probs=pd.DataFrame(tempubox,columns=["gene","pos","probability"])
    lbox_probs=lbox_probs.drop(0)
    ubox_probs=ubox_probs.drop(0)

    lbox_probs.to_csv(outPrefix.strip()+"_lrnafold.txt",sep="\t")
    ubox_probs.to_csv(outPrefix.strip()+"_urnafold.txt",sep="\t")
def main():
    parser=creat_parse()
    (options,args)=parser.parse_args()
    (input,outPrefix,fastaFile)=(options.input,options.outPrefix,options.fastaFile)
    print("Start preprocessing the data...",file=sys.stderr)
    preprocess(input,outPrefix,fastaFile)
    print("Finish the step of processing...")


if __name__=="__main__":
    main()

## 脚本使用
python ProcessingRNAfold.py -i inputDir -o outprefix -f trans.fa

## 其中inputDir为RNAfold输出结果，ouprefix为前缀，trans.fa为用户感兴趣的每个序列
## 输出结果包含两个文件： outprefix_down_lrnafold.txt,outprefix_rrnafold.txt

```

+ 绘制二级结构可能性的metagene plot

```
#!/usr/bin/python
# -*- coding: UTF-8 -*-

import pandas as pd
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


inputFile=sys.argv[1]
left_pos=sys.argv[2]
right_pos=sys.argv[3]


def extractSpecificRangeScores(inputFile,left_pos,right_pos):
    lbox_data=pd.read_table(inputFile,sep='\t')
    geneList=np.unique(lbox_data.iloc[:,1])
    prob={}
    for gene in geneList:
        # tmp=lbox_data[lbox_data['gene']==gene]
        tmp=lbox_data.loc[lbox_data.iloc[:,1]==gene,:]
        tmp=tmp.iloc[np.arange(int(left_pos)-1,int(right_pos)),:]
        tmp.iloc[:,3]=np.where(tmp.iloc[:,3]>=0.5,1,0)
        prob[gene]=tmp.iloc[:,3].values
    prob=pd.DataFrame(prob)
    final_prob=prob.apply(sum,axis=1)
    return prob,final_prob
def plot_final_score(final_prob):
    x_range=np.arange(len(final_prob))
    fig=plt.figure()
    plt.plot(x_range,final_prob.values)
    plt.xlabel("Position (nt)")
    plt.ylabel("Probability scores")
    plt.show()


def main():
    prob,final_prob=extractSpecificRangeScores(inputFile,left_pos,right_pos)
    plot_final_score(final_prob)



if __name__ =="__main__":
    main()

```