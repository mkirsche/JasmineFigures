BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
SURVIVOR_PATH='/home/mkirsche/git/SURVIVOR'
JASMINEEVALSRCPATH='/home/mkirsche/git/JasmineResults/src'
OUTPREFIX='population'
FILELIST=$BINDIR/$OUTPREFIX.filelist.txt
#$SURVIVOR_PATH/Debug/SURVIVOR merge $FILELIST 1000 1 1 1 0 1 $BINDIR/$OUTPREFIX.survivor.vcf

javac $JASMINEEVALSRCPATH/ShuffleFile.java
javac $JASMINEEVALSRCPATH/VcfToTsv.java

for i in {1..10}
do
 echo $i
 if [ -d $BINDIR/shuffle$i ]
 then
   rm -r $BINDIR/shuffle$i
 fi
 mkdir $BINDIR/shuffle$i

 NEWFILELIST=$BINDIR/shuffle$i/$OUTPREFIX.filelist.txt
 
 java -cp $JASMINEEVALSRCPATH ShuffleFile in_file=$FILELIST out_file=$NEWFILELIST seed=$i
 
 mergedvcf=$BINDIR/shuffle$i/$OUTPREFIX.survivor.vcf
 $SURVIVOR_PATH/Debug/SURVIVOR merge $BINDIR/shuffle$i/$OUTPREFIX.filelist.txt 1000 1 1 1 0 1 $mergedvcf
 mergedtsv=$BINDIR/shuffle$i/$OUTPREFIX.survivor.tsv
 java -cp $JASMINEEVALSRCPATH VcfToTsv $mergedvcf $mergedtsv
done

