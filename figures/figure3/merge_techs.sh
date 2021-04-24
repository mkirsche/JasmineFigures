BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
MAX_DIST=50
ALIGNER='winnowmap'
SRCPATH=$BINDIR'/../../src'
JASMINEPATH=$BINDIR'/../../Jasmine'

PBFILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/pb/'$ALIGNER'/HG002v'*'_'$MAX_DIST'md_'*.ism.vcf`
ONTFILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/ont/'$ALIGNER'/'*'_'$MAX_DIST'md_'*.ism.vcf`
HIFIFILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/hifi/'$ALIGNER'/'*'_'$MAX_DIST'md_'*.ism.vcf`

echo 'PB file: '$PBFILE
echo 'ONT file: '$ONTFILE
echo 'Hifi file: '$HIFIFILE

prefix='hg002_crosstech'
filelist=$BINDIR/$prefix.filelist.txt
echo $HIFIFILE > $filelist
echo $PBFILE >> $filelist
echo $ONTFILE >> $filelist

echo 'File list: '$filelist

mergedvcf=$BINDIR/$prefix.merged.vcf

$JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf max_dist_linear=0.5 min_dist=100
$JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
echo 'Merged VCF: '$mergedvcf

mergedtsv=$BINDIR/$prefix.merged.tsv
java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv
echo 'Merged TSV: '$mergedtsv
