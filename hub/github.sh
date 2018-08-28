git add -A
git commit -a
git push

########
ddir=~/workspace/9.NT-ChIP/2.public/3.bowtie2/c.merge_rep
wdir=~/git/discovery/hub/2.HM
for i in $ddir/[mIT]*bam; do 
nohup bamCoverage -b $i -o $wdir/$(basename $i .sorted.bam).bw --normalizeUsing RPKM -p 8 --ignoreDuplicates --smoothLength 40000 -bs 20000 \
--samFlagInclude 2 --minFragmentLength 100 --maxFragmentLength 500 &
done



http://genome-asia.ucsc.edu/cgi-bin/hgTracks?udcTimeout=1&db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/1.ES/hub.txt
http://genome-asia.ucsc.edu/cgi-bin/hgTracks?udcTimeout=1&db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/3.compartment/hub.txt
http://genome-asia.ucsc.edu/cgi-bin/hgTracks?udcTimeout=1&db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/2.HM/hub.txt

http://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/1.ES/hub.txt
http://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/3.compartment/hub.txt
http://genome-asia.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/2.HM/hub.txt



for k in zygote 2cell 4cell 8cell Morula ICM TE ; do 
sed 's|8cell|'$k'|' track.temple >> 2.HM/track.txt;done


############################################Insulation score
wdir=~/git/discovery/hub/4.IS
ddir=~/workspace/8.NT-HiC/f.IS_ALL/2.except_res40k_is1M_ids200k_nt01/4.bw
mkdir -p $wdir
cp $ddir/*bw $wdir

cp $wdir/../3.compartment/*txt $wdir

#creat track.temple
rm $wdir/track.txt
for k in ${keys[@]}; do 
sed 's|CC|'$k'|' $wdir/track.temple >> $wdir/track.txt;
done

http://genome-asia.ucsc.edu/cgi-bin/hgTracks?udcTimeout=1&db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/4.IS/hub.txt

############################################Directional index
wdir=~/git/discovery/hub/5.DI
ddir=~/workspace/8.NT-HiC/g.DI_ALL/2.except_norm_100M_20180825/4.directional_index
mkdir -p $wdir

for i in ${keys[@]}; do 
bedSort $ddir/$i.bedGraph $wdir/$i.bedGraph
bedGraphToBigWig $wdir/$i.bedGraph ~/ann/mm10.chromSizes $wdir/$i.bw
rm $wdir/$i.bedGraph
done


cp $wdir/../4.IS/*txt $wdir

#creat temple
rm $wdir/track.txt
for k in ${keys[@]}; do 
sed 's|CC|'$k'|' $wdir/temple.txt >> $wdir/track.txt;
done


http://genome-asia.ucsc.edu/cgi-bin/hgTracks?udcTimeout=1&db=mm10&hubUrl=https://raw.githubusercontent.com/rysterzhu/discovery/master/hub/4.IS/hub.txt
