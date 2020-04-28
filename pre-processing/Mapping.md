##**Mapping**

目前用于序列比对的工具有很多包括BWA, bowtie, bowtie2, tophat, hisat2, STAR等等，每种方法各有优缺点，这里不再详细讨论，具体可参见各自工具的说明文档，我们平时用bowtie/2和STAR比较多，所以这里会以这两个工具为例，介绍其简单的用法，详细的参数见其说明文档，这里也不再说明。

###**bowtie使用**

第一步： 安装[bowtie](http://bioconda.github.io/recipes/bowtie/README.html)

```
## conda 下载
conda install bowtie

## 其他方式见官网
```

第二步： 建立索引

```
## human_genome为索引的前缀名
fastaGenome=/workdata/LabMember2/xdxing/GenomeFiles/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
bowtie-build -f $fastaGenome human_genome

```
第三步： mapping

```
## 去除rRNA污染
bowtie_noncoding_index=/workdata/LabMember2/lifj/lifj/data/Reference/human/rRNA_bowtieIndex/human_rRNA
fastqFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE63591/03.filter
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE63591/05.contam

bowtie -n 0 -y -a --norc --best --strata -S -p 4 -l 15 --un=$workdir/noncontam_$i.fastq $bowtie_noncoding_index -q $fastqFile/$i.trimmed.Qfilter.fastq $workdir/$i.alin
```

具体的参数细节详见[bowtie官网](http://bowtie-bio.sourceforge.net/manual.shtml)

###**bowtie2使用**

这里不再详述，具体细节见官网[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

##**STAR使用**

第一步：安装[STAR](https://github.com/alexdobin/STAR)

```
## download 

# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
cd STAR-2.7.3a

# Alternatively, get STAR source using git
git clone https://github.com/alexdobin/STAR.git


## Compile
cd STAR/source
make STAR

```

第二步： 建立索引,详见[文档](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

```
fastaGenome=/workdata/LabMember2/xdxing/GenomeFiles/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=/workdata/LabMember2/xdxing/GenomeFiles/Homo_sapiens/Ensembl/GRCh38/Homo_sapiens.GRCh38.88.gtf
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./STAR_Human_Ensembl_GRCh38_Ensembl --genomeFastaFiles $fastaGenome --sjdbGTFfile $gtf
```

第三步： STAR mapping,详见[文档](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

```
STAR_genome_index=/workdata/LabMember2/lifj/lifj/data/Reference/human/STAR_Human_Ensembl_GRCh38_Ensembl
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE63591/07.STAR
fastqFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE63591/05.contam


for i in SRR16614{91..98} SRR20754{17..20};
do
    mkdir -p $workdir/${i}_STAR
    
    ##  mapping for single end (SE)
    STAR --runThreadN 8 --outFilterType Normal --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir $STAR_genome_index --readFilesIn $fastqFile/noncontam_$i.fastq --outFileNamePrefix  $workdir/${i}_STAR/$i. --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All
    
    ## mapping for pair end (PE)
    STAR --runThreadN 8 --outFilterType Normal --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapNmax 1 --genomeDir $STAR_genome_index --readFilesIn read_1.fa read_2.fa --outFileNamePrefix  $workdir/${i}_STAR/$i. --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All
done

```

第四步： sort bam files and index

```
## sort
for i in SRR16614{91..98} SRR20754{17..20};
do
    samtools sort -T $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted -o $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.bam
done

## set index
for i in SRR16614{91..98} SRR20754{17..20};
do
    samtools index  $workdir/${i}_STAR/$i.Aligned.toTranscriptome.out.sorted.bam
    samtools index $workdir/${i}_STAR/$i.Aligned.sortedByCoord.out.bam
done
```

第五步（optional）: 构建bigwig文件，方便IGV查看

```
workdir=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE70211/07.STAR
bamFile=/workdata/LabMember2/lifj/lifj/Project/05.Ribo_seq_human/GSE70211/07.STAR

for i in SRR20759{25..34} SRR20759{36..39}
do
## RNA
bamCoverage -b $bamFile/${i}_STAR/$i.Aligned.sortedByCoord.out.bam -o $workdir/${i}_STAR/${i}.bw --normalizeUsing CPM --binSize 1
## RPF

done

```