BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
OUTPREFIX='population'
FILELIST=$BINDIR'/'$OUTPREFIX'.filelist.txt'
SRCPATH=$BINDIR'/../../src'

SHUFFLEDIR=$BINDIR/shuffle
if [ ! -d $SHUFFLEDIR ]
then
  mkdir $SHUFFLEDIR
fi

cd $SHUFFLEDIR

for i in {1..1}
do
 echo $i
 NEWFILELIST=$SHUFFLEDIR/$OUTPREFIX'_'$i.filelist.txt
 
 java -cp $SRCPATH ShuffleFile in_file=$FILELIST out_file=$NEWFILELIST seed=$i
 
 # Remove extension from filelist for making caller-specific versions of it
 filelist_noext=`echo "${NEWFILELIST%.*}"`

 fixedlist=$filelist_noext'.fixed.txt'
 if [ -r $fixedlist ]
 then
   rm $fixedlist
 fi

 javac $BINDIR/../figure2/*.java

 echo 'Adding strand bias headers'
 for filename in `cat $NEWFILELIST`
 do
   noext=`echo "${filename%.*}"`
   outfile=$noext.fixed.vcf 
   echo $filename $outfile $fixedlist
   echo $outfile >> $fixedlist
   java -cp $BINDIR/../figure2 AddStrandBiasHeader $filename $outfile
 done
 echo 'Done adding strand bias headers'
 
 prefix='population_shuffle_'$i
 discsuppvec='0'
 $BINDIR/../figure2/SvPopulationAnalysis/runmerging.sh $fixedlist $prefix $discsuppvec
done
