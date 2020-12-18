if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

filelist=$BINDIR/jasmine_md50.filelist.txt
prefix='hg002_hifi'
discsuppvec='100'

# Remove extension from filelist for making caller-specific versions of it
filelist_noext=`echo "${filelist%.*}"`

fixedlist=$filelist_noext'.fixed.txt'
if [ -r $fixedlist ]
then
  rm $fixedlist
fi

javac $BINDIR/*.java

for filename in `cat $filelist`
do
  noext=`echo "${filename%.*}"`
  outfile=$noext.fixed.vcf 
  echo $outfile >> $fixedlist
  java -cp $BINDIR AddStrandBiasHeader $filename $outfile
done

$BINDIR/SvPopulationAnalysis/runmerging.sh $fixedlist $prefix $discsuppvec




