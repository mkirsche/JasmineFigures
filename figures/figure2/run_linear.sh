if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi
JASMINEPATH=$BINDIR'/../../Jasmine'
SRCPATH=$BINDIR'/../../src'
filelist=$BINDIR/crosstool.filelist.txt
sonfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg002/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
fatherfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg003/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
motherfile=`ls $BINDIR/../../schatz_pipeline/ash_trio/hg004/hifi/winnowmap/*vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf`
echo $sonfile > $filelist
echo $fatherfile >> $filelist
echo $motherfile >> $filelist
SLOPE=$1
MIN=$2
MAX=2000
SLOPE_DECIMAL=$(printf %.3f $(echo "$SLOPE/100" | bc -l));
prefix='jasmine_linear'$SLOPE'_'$MIN
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

mergedvcf=$BINDIR/$prefix.merged.vcf
$JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf max_dist_linear=$SLOPE_DECIMAL min_dist=$MIN --mutual_distance
$JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf

echo 'Merged VCF: '$mergedvcf

mergedtsv=$BINDIR/$prefix.merged.tsv
javac $SRCPATH/*.java
java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv
echo 'Merged TSV: '$mergedtsv



