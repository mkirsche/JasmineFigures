if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

filelist=$BINDIR/crosstool.filelist.txt
sonfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg002/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
fatherfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg003/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
motherfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg004/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
echo $sonfile > $filelist
echo $fatherfile >> $filelist
echo $motherfile >> $filelist
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

$BINDIR/SvPopulationAnalysis/run_svpop/run_svpop.sh $fixedlist $prefix $discsuppvec


