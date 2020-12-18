BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
SURVIVOR_PATH='/home/mkirsche/git/SURVIVOR'
JASMINEEVALSRCPATH='/home/mkirsche/git/JasmineResults/src'
TABLESRCPATH='/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/src'
OUTPREFIX='population'
FILELIST=$BINDIR/$OUTPREFIX.filelist.txt

javac $JASMINEEVALSRCPATH/ShuffleFile.java
javac $JASMINEEVALSRCPATH/VcfToTsv.java

tablefile=$BINDIR/$OUTPREFIX.counts.txt

echo ITERATION$'\t'COUNT$'\t'MERGER > $tablefile

SHUFFLEDIR=$BINDIR/shuffle
if [ ! -d $SHUFFLEDIR ]
then
  mkdir $SHUFFLEDIR
fi

for i in {1..10}
do
 echo $i

 NEWFILELIST=$SHUFFLEDIR/$OUTPREFIX'_'$i.filelist.txt
 
 java -cp $JASMINEEVALSRCPATH ShuffleFile in_file=$FILELIST out_file=$NEWFILELIST seed=$i
 
 survivorvcf=$SHUFFLEDIR/$OUTPREFIX'_'$i.survivor.vcf
 $SURVIVOR_PATH/Debug/SURVIVOR merge $NEWFILELIST 1000 1 1 1 0 1 $survivorvcf
 
 # Build merging table from SURVIVOR
 echo 'Building SURVIVOR merging table'
 survivorsmalltsv=$SHUFFLEDIR/$OUTPREFIX'_'$i.survivor_small.tsv
 java -cp /home/mkirsche/eclipse-workspace/Jasmine/src:$TABLESRCPATH BuildMergingTable vcf_file=$survivorvcf out_file=$survivorsmalltsv vcf_filelist=$NEWFILELIST mode=survivor
 
 # Add INFO annotations to SURVIVOR merging table
 echo 'Annotating SURVIVOR merging table'
 survivortsv=$SHUFFLEDIR/$OUTPREFIX'_'$i.survivor.tsv
 java -cp /home/mkirsche/eclipse-workspace/Jasmine/src:$TABLESRCPATH AugmentMergingTable table_file=$survivorsmalltsv out_file=$survivortsv vcf_filelist=$NEWFILELIST
 
 survivorcount=`cat $survivortsv | awk -v cols='SPECIFIC_FLAG,PRECISE_FLAG' 'BEGIN {
   FS=OFS="\t"
   nc=split(cols, a, ",")
}
NR==1 {
   for (i=1; i<=NF; i++)
      hdr[$i]=i
}
{
   for (i=1; i<=nc; i++)
      if (a[i] in hdr)
         printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)
}' | awk '{ if(($1 == 1) && ($2 == 1)) { print } }' | wc -l`
 
 
 rm $survivorvcf
 rm $survivorsmalltsv
 rm $survivortsv
 
  echo 'SURVIVOR produced '$survivorcount' variants'
  echo $i$'\t'$survivorcount$'\t'SURVIVOR >> $tablefile
 
 # Run Jasmine
 echo 'Running Jasmine'
 jasminevcf=$SHUFFLEDIR/$OUTPREFIX'_'$i.jasmine.vcf
 ~/eclipse-workspace/Jasmine/jasmine file_list=$NEWFILELIST out_file=$jasminevcf

 # Build merging table from Jasmine
 echo 'Building Jasmine merging table'
 jasminesmalltsv=$SHUFFLEDIR/$OUTPREFIX'_'$i.jasmine_small.tsv
 java -cp /home/mkirsche/eclipse-workspace/Jasmine/src:$TABLESRCPATH BuildMergingTable vcf_file=$jasminevcf out_file=$jasminesmalltsv vcf_filelist=$NEWFILELIST
 
 # Add INFO annotations to Jasmine merging table
 echo 'Annotating Jasmine merging table'
 jasminetsv=$SHUFFLEDIR/$OUTPREFIX'_'$i.jasmine.tsv
 java -cp /home/mkirsche/eclipse-workspace/Jasmine/src:$TABLESRCPATH AugmentMergingTable table_file=$jasminesmalltsv out_file=$jasminetsv vcf_filelist=$NEWFILELIST
 
 jasminecount=`cat $jasminetsv | awk -v cols='SPECIFIC_FLAG,PRECISE_FLAG' 'BEGIN {
   FS=OFS="\t"
   nc=split(cols, a, ",")
}
NR==1 {
   for (i=1; i<=NF; i++)
      hdr[$i]=i
}
{
   for (i=1; i<=nc; i++)
      if (a[i] in hdr)
         printf "%s%s", $hdr[a[i]], (i<nc?OFS:ORS)
}' | awk '{ if(($1 == 1) && ($2 == 1)) { print } }' | wc -l`
 
 rm $jasminevcf
 rm $jasminesmalltsv
 rm $jasminetsv
 
 echo 'Jasmine produced '$jasminecount' variants'
 echo $i$'\t'$jasminecount$'\t'Jasmine >> $tablefile

done

