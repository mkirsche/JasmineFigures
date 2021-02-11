BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
ALIGNER='winnowmap'
md='50'
JASMINEPATH=$BINDIR'/../../../Jasmine'
SRCPATH=$BINDIR'/../../../src'


dists=( 100 250 500 750 1000 1500 2000 3000 4000 5000 )
for dist in ${dists[@]}
do
    echo $dist
    HG002FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg002/hifi/'$ALIGNER'/'*.ism.vcf`
    HG003FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg003/hifi/'$ALIGNER'/'*.ism.vcf`
    HG004FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg004/hifi/'$ALIGNER'/'*.ism.vcf`
    echo 'HG002 variants: '$HG002FILE
    echo 'HG003 variants: '$HG003FILE
    echo 'HG004 variants: '$HG004FILE
    
    prefix='jasmine_dist'$dist
    filelist=$BINDIR/$prefix.filelist.txt
    echo $HG002FILE > $filelist
    echo $HG003FILE >> $filelist
    echo $HG004FILE >> $filelist

    echo 'File list: '$filelist
    
    mergedvcf=$BINDIR/$prefix.merged.vcf

    $JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf max_dist=$dist
    $JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
    
    echo 'Merged VCF: '$mergedvcf
    
    mergedtsv=$BINDIR/$prefix.merged.tsv
    javac $SRCPATH/*.java
    java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv
    echo 'Merged TSV: '$mergedtsv
done

Rscript $BINDIR/../../jasmine_plot/jasmine_dist_discordance.R
