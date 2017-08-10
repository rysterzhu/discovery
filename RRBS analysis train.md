## <div align="center"> RRBS analysis train</div>

1. 使用/home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/的T7和pT7数据练习；
2. 使用trim_galore + bismark + methylKit pipline

![](http://www.bioinformatics.babraham.ac.uk/images/babraham_bioinformatics.gif)
* [babraham Methylation_Course](https://www.bioinformatics.babraham.ac.uk/training/Methylation_Course/)

---
#### 1.数据准备

	 3140  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.1/5-2_L8_7046.R2.clean.fastq.gz p1.R2.fastq.gz
	 3141  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.2/5-3_L8_7047.R1.clean.fastq.gz p2.R1.fastq.gz
	 3142  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.2/5-3_L8_7047.R2.clean.fastq.gz p2.R2.fastq.gz
	 3143  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.3/5-1_L2_7045.R1.clean.fastq.gz p31.R1.fastq.gz
	 3144  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.3/5-1_L2_7045.R2.clean.fastq.gz p31.R2.fastq.gz
	 3145  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.3/5-1_L3_A031.R1.clean.fastq.gz p32.R1.fastq.gz
	 3146  ln -s /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/pT7OC.3/5-1_L3_A031.R2.clean.fastq.gz p32.R2.fastq.gz
	 3150  for i in /home/share/DATA/clean_reads/GaoSR_mouse_oocyte_RRBS-Seq_20170626/T7OC.*/*;do  j=${i%/*};j=${j##*.}; k=${i##*.R};k=${k%%.*}; ln -s $i t$j.R$k.fastq.gz; done
![](http://i.imgur.com/G3kp9rl.jpg)
---
#### 2. fastqc
	[~/workspace/5.rrbs_study/0.data]$nohup fastqc -o 1.qc/ 0.links/* > 1.qc/run.log &
	[~/workspace/5.rrbs_study/0.data/1.qc]$bash ~/codes/qc_summary.sh
	awk -v FS="\t" -v OFS=" | " '{$1=$1;print } NR==1{for(i=1;i<=NF;i++){$i="--"};print}' states_1x.tab # 将table转markdown表格格式

> T7.3 R1质量不好，其他都还行
> 但仍然需要所有的样本都进行trim，因为RRBS两端可能有CT偏差


sample | total | length | %GC | Q20 | Q30 
--- | --- | --- | --- | --- | ---
p1.R2 | 9760244 | 125 | 32 | 0.993167179017246 | 0.953253217849882
p2.R1 | 11867309 | 125 | 34 | 0.996014429218958 | 0.942492185886455 
p2.R2 | 11867309 | 125 | 33 | 0.99257312672991 | 0.951581778143638
p31.R1 | 5462552 | 125 | 35 | 0.986030705062396 | 0.851489010997058 
p31.R2 | 5462552 | 125 | 34 | 0.990368055077554 | 0.915194949173939
p32.R1 | 6369957 | 125 | 35 | 0.993463849128024 | 0.917539003795473 
p32.R2 | 6369957 | 125 | 34 | 0.989659427842292 | 0.916141349148825
t1.R1 | 9880816 | 125 | 34 | 0.987846044294317 | 0.853355431373279
t1.R2 | 9880816 | 125 | 33 | 0.990461516538715 | 0.922468549156264
t2.R1 | 7273755 | 125 | 34 | 0.981606886676826 | 0.7844187768216 
t2.R2 | 7273755 | 125 | 33 | 0.987218018753725 | 0.903860798170958
t3.R1 | 4214607 | 125 | 35 | 0.963983830520853 | 0.659060738047462 
t3.R2 | 4214607 | 125 | 34 | 0.976846239756162 | 0.830806525970274

Sample | Basic Statistics | Per base sequence quality | Per tile sequence quality | Per sequence quality scores | Per base sequence content | Per sequence GC content | Per base N content | Sequence Length Distribution | Sequence Duplication Levels | Overrepresented sequences | Adapter Content | Kmer Content
-- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | --
p1.R2_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
p2.R1_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
p2.R2_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
p31.R1_fastqc | PASS | PASS | WARN | PASS | FAIL | WARN | PASS | PASS | FAIL | FAIL | PASS | FAIL
p31.R2_fastqc | PASS | PASS | FAIL | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
p32.R1_fastqc | PASS | PASS | WARN | PASS | FAIL | WARN | PASS | PASS | FAIL | FAIL | PASS | FAIL
p32.R2_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | WARN
t1.R1_fastqc | PASS | WARN | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
t1.R2_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL
t2.R1_fastqc | PASS | PASS | FAIL | PASS | FAIL | WARN | PASS | PASS | FAIL | FAIL | PASS | FAIL
t2.R2_fastqc | PASS | PASS | FAIL | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | WARN
t3.R1_fastqc | PASS | FAIL | FAIL | PASS | FAIL | WARN | PASS | PASS | FAIL | FAIL | PASS | WARN
t3.R2_fastqc | PASS | PASS | WARN | PASS | FAIL | FAIL | PASS | PASS | FAIL | FAIL | PASS | FAIL

![](http://i.imgur.com/U312Ftz.jpg)
---
#### 3.trim
	nohup trim_galore --illumina --fastqc --rrbs --paired -o 2.trim/ 0.links/*fastq.gz > 2.trim/trim.log &
	[~/workspace/5.rrbs_study/0.data]$for i in 0.links/*R1.fastq.gz; do
	> nohup trim_galore --illumina --fastqc --rrbs --paired --length 50 -o 2.trim/ $i ${i/R1/R2} >> 2.trim/trim.log &
	> done
> 默认控制的质量是20(-q 20)，去掉的reads不多
> [trim_galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/):
> remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing)
> 集成了cutadapt 和 fastqc

sample | total | length | %GC | Q20 | Q30
--- | --- | --- | --- | --- | --- 
p1.R1_val_1 | 9760181 | 50-125 | 33 | 1 | 0.964745530846201 
p1.R2_val_2 | 9760181 | 50-123 | 32 | 0.999999795085767 | 0.963135007434801
p2.R1_val_1 | 11867272 | 50-125 | 34 | 0.999999831469299 | 0.96660816966053
p2.R2_val_2 | 11867272 | 50-123 | 33 | 0.999999831469285 | 0.961808327827294
p31.R1_val_1 | 5462432 | 50-125 | 34 | 0.999999633862853 | 0.912728465002964 
p31.R2_val_2 | 5462432 | 50-123 | 34 | 0.999999633862853 | 0.933397822289478
p32.R1_val_1 | 6369850 | 51-125 | 35 | 0.999999686020962 | 0.940591555345399
p32.R2_val_2 | 6369850 | 50-123 | 34 | 0.999999686020913 | 0.931278141356839
t1.R1_val_1 | 9880691 | 50-125 | 34 | 0.999999797585048 | 0.927904753239474 
t1.R2_val_2 | 9880691 | 50-123 | 32 | 0.99999969637751 | 0.941334062567082
t2.R1_val_1 | 7262971 | 50-125 | 34 | 0.999999311576583 | 0.825811114238083 
t2.R2_val_2 | 7262971 | 50-123 | 33 | 0.999993391134289 | 0.921704905609564
t3.R1_val_1 | 4212408 | 50-125 | 35 | 0.999995014751587 | 0.725989619288448 
t3.R2_val_2 | 4212408 | 50-123 | 33 | 0.999997388667005 | 0.860081691991849

![](http://i.imgur.com/YrCItDG.jpg)
---
#### 4.bismark 
* [bismark 原理](http://blog.csdn.net/qq_29300341/article/details/53302961)
* [Bismark_User_Guide](https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html "Bismark_User_Guide")
##### 1)准备bismark需要转换的index
	bismark_genome_preparation --verbose /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/
	bismark_genome_preparation --verbose /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/ &
> bismark需要使用bowtie进行比对，默认创建的bowtie2的index；使用bowtie1需要--bowtie1指定bowtie位置
得到CT_conversion/  GA_conversion/两种转换的indexes
* 已统一将bismark indexes建立链接到/home/share/bismark_index
##### 2)默认的bismark
	for i in ~/workspace/5.rrbs_study/0.data/2.trim/*R1_val_1.fq.gz; do
	o=${i##*/};o=${o%%.*}
	nohup bismark -p 8 -N 1 /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ -1 $i -2 ${i/R1_val_1/R2_val_2} -o $o > logs/$o.log &
	done 
> bismark有bowtie2的大部分参数
默认使用的bowtie2参数：-q -N 1 --score-min L,0,-0.2 -p 8 --reorder --ignore -quals --no-mixed --no-discordant --dovetail --maxins 500 --norc
默认directional library（只比对OT和OB链，不比对CTOT和CTOB链）
-B basename

##### 3)测试directional和 non-directional library


	nohup bismark -N 1 /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ -1 ~/workspace/5.rrbs_study/0.data/t11.R1.fastq.gz -2 ~/workspace/5.rrbs_study/0.data/t11.R2.fastq.gz -o t11.bismark1 > t11.bismark1.log &
	nohup bismark -N 1 /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ -1 ~/workspace/5.rrbs_study/0.data/t11.R1.fastq.gz -2 ~/workspace/5.rrbs_study/0.data/t11.R2.fastq.gz -o t11.bismark3 --non_directional > t11.bismark3.log &
	nohup bismark_methylation_extractor --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ t11.R1_bismark_bt2_pe.bam -o extractor2/ >> methylation_extractor.log &

> 有--non_directional时：CTOT和CTOB的比对率只有0.02%；
> CTOT和CTOB的CpG只有200多个
> 判断应该是directional的建库
具体差异细节见两个比对的日志 t11.bismark1.log t11.bismark3.log

---
#### 5.reports
	bismark2summary *bam
	bismark2report --alignment_report p1.R1_val_1_bismark_bt2_PE_report.txt
> bismark2report会使用deduplication report file（RRBS不用考虑）、splitting report file、M-bias report file、nucleotide coverage report file，可用参数指定，不指定将寻找当前目录下所有符合的report
![](http://i.imgur.com/8PxQDKb.jpg)
[p1.all.report](p1.R1_val_1_bismark_bt2_PE_report.html)
[bismark_summary_report](bismark_summary_report.html)

##### -总体情况
![](http://i.imgur.com/RzFE3xi.jpg)
![](http://i.imgur.com/EaWXAqD.jpg)

##### -T7.1的甲基化
![](http://i.imgur.com/Vxy8LUO.jpg)
> 比对率（T7.1）
![](http://i.imgur.com/G3NjOCz.jpg)
> 甲基化的占比CpG（T7.1）

##### -bam2nuc （统计(di-)nucleotide）
	bam2nuc --genome_folder /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ ../p1.R1_val_1_bismark_bt2_pe.bam

(di-)nucleotide | count sample | percent sample | count genomic | percent genomic | coverage
-- | -- | -- | -- | -- | --
A | 232548822 | 23.01 | 744681828 | 29.11 | 0.312
C | 250105490 | 24.75 | 534146040 | 20.88 | 0.468
G | 307300011 | 30.41 | 534300392 | 20.88 | 0.575
T | 220647722 | 21.83 | 745397519 | 29.13 | 0.296
AA | 53445908 | 5.33 | 232348065 | 9.08 | 0.230
AC | 47473139 | 4.74 | 136264156 | 5.33 | 0.348
AG | 90045519 | 8.98 | 188025079 | 7.35 | 0.479
AT | 41281559 | 4.12 | 188044390 | 7.35 | 0.220
CA | 69161868 | 6.90 | 190623924 | 7.45 | 0.363
CC | 67525762 | 6.74 | 134073230 | 5.24 | 0.504
CG | 35656813 | 3.56 | 21342779 | 0.83 | 1.671
CT | 72094995 | 7.19 | 188105878 | 7.35 | 0.383
GA | 73883127 | 7.37 | 159020818 | 6.22 | 0.465
GC | 75146249 | 7.50 | 104696885 | 4.09 | 0.718
GG | 99256956 | 9.90 | 134087829 | 5.24 | 0.740
GT | 57685988 | 5.76 | 136494721 | 5.33 | 0.423
TA | 35013312 | 3.49 | 162688871 | 6.36 | 0.215
TC | 54968450 | 5.48 | 159111635 | 6.22 | 0.345
TG | 80700184 | 8.05 | 190844482 | 7.46 | 0.423
TT | 48974415 | 4.89 | 232752388 | 9.10 | 0.210
---
#### 6.bismark_methylation_extractor
	nohup bismark_methylation_extractor --gzip --bedGraph --buffer_size 10G --cytosine_report --genome_folder /home/share/Illumina_iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/ *bam -o extractor/ >> methylation_extractor.log &
> 得到CpG CHG CHH的甲基化位点，所有甲基化位点(默认是CpG的，可指定)的coverage、ratio，一些报告
--ignore 2(--ignore_r2 2 )：切除reads 5‘两个碱基
--comprehensive：合并四种链
--bedGraph：输出methy ratio
--cytosine_report --genome_folder：输出所有CpG的甲基化（.CpG_report.txt.gz），需指定genome
--CX：指定CpG CHH CHG
![](http://i.imgur.com/V5uA2Rr.jpg)
> 3:methy ratio; 4:methyC; 5: Cytosine;

---
#### 7.methylKit
* [methylKit user guide](http://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html#2_basics)  
##### 1)read files
```R
file.list=as.list(list.files("~/workspace/5.rrbs_study/1.align/t-p/extractor",pattern = "*cov.gz",full.names = T))

met_raw=methRead(file.list,
           sample.id=list("treat1","treat2","treat3","treat4","control1","control2","control3"),
           assembly="mm9",
           treatment=c(1,1,1,1,0,0,0),
           context="CpG",
          # dbtype = "tabix",
          # dbdir = "methylDB",
           pipeline = "bismarkCoverage", #("amp", "bismark","bismarkCoverage", "bismarkCytosineReport" or a list)
          header = FALSE,
          mincov = 3 #(minimum read coverage,defaults to 10)
           )

met_raw2 <- normalizeCoverage(met_raw)
meth=unite(met_raw, destrand=FALSE)
meth2=unite(met_raw2,destrand=FALSE,
            min.per.group=3L #(minimum number of samples per replicate needed to cover a region/base,default all)
)
```
##### 2) plots
	getMethylationStats(met_raw[[1]],plot=TRUE,both.strands=FALSE)
![](http://i.imgur.com/R98lUsK.jpg)

	getCoverageStats(met_raw[[1]],plot=TRUE,both.strands=FALSE)
![](http://i.imgur.com/YsQgaqT.jpg)

	clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
![](http://i.imgur.com/xiZ1r2R.jpg)

	getCorrelation(meth,method="pearson",plot=T)
![](http://i.imgur.com/D6Tela8.jpg)

	PCASamples(meth)
![](http://i.imgur.com/90vLtLf.jpg)
> normalized
##### 3) calculate Differentially Methylation
```r
myDiff=calculateDiffMeth(meth)
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.05)

diffMethPerChr(myDiff,plot=T,qvalue.cutoff=0.01, meth.cutoff=25)
```

![](http://i.imgur.com/9xHYliP.png)
![](http://i.imgur.com/XuR5q4r.jpg)
##### 4) Annotation
```r
library(genomation)
gene.obj=readTranscriptFeatures("/home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.RefSeq.bed")
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank("/home/share/ann/useful_annotation/download_FTP/UCSC_table_broswer/mm9.CGI.bed",
                           feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")


promoters=regionCounts(met_raw,gene.obj$promoters)

head(promoters[[1]])

diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)
plotTargetAnnotation(diffAnn,precedence=TRUE,
    main="differential methylation annotation")
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="differential methylation annotation")
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)
```
![](http://i.imgur.com/iG3gOsM.jpg)

	plotTargetAnnotation(diffAnn,precedence=TRUE,main="differential methylation annotation")
![](http://i.imgur.com/Vru8Om0.jpg)

	plotTargetAnnotation(diffCpGann,col=c("green","gray","white"), main="differential methylation annotation")
![](http://i.imgur.com/E9O6D91.jpg)


----
#### 8.lambda测试
	nohup bismark -p 8 -N 1 /home/share/bismark_index/lambda -1 ~/workspace/5.rrbs_study/0.data/2.trim/p1.R1_val_1.fq.gz -2 ~/workspace/5.rrbs_study/0.data/2.trim/p1.R2_val_2.fq.gz -o p1 > p1.log &

![](http://i.imgur.com/ftM2v6v.jpg)
> lambda转化率 = (125244+137007+143947)/409667 = 99.1%

	for i in ~/workspace/5.rrbs_study/0.data/2.trim/*R1_val_1.fq.gz; do
	o=${i##*/};o=${o%%.*}
	nohup bismark -p 8 -N 1 /home/share/bismark_index/lambda -1 $i -2 ${i/R1_val_1/R2_val_2} -o $o > logs/$o.log &
	done 
	for i in *; do echo ${i%%.*} `awk -F ":\t" '$0~/^Total number/{a=$2} $0~/^Total unmethylated/{b+=$2} END{print b/a}' $i`; done

sample| convert ratio
--|--
p1| 0.991532
p2| 0.994595
p31| 0.993302
p32| 0.993291
t1| 0.994703
t2| 0.994134
t3| 0.993318



----
