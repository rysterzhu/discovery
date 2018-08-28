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


