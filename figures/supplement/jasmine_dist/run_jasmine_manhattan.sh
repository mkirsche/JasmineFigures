BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
ALIGNER='winnowmap'
md='50'
JASMINEPATH=$BINDIR'/../../../Jasmine'
SRCPATH=$BINDIR'/../../../src'

dists=( 1000 )
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
    
    manhattanvcf=$BINDIR/'jasmine_manhattan'.merged.vcf
    euclideanvcf=$BINDIR/'jasmine_euclidean'.merged.vcf

    $JASMINEPATH/jasmine file_list=$filelist out_file=$manhattanvcf max_dist=$dist kd_tree_norm=1
    $JASMINEPATH/jasmine file_list=$filelist out_file=$euclideanvcf max_dist=$dist kd_tree_norm=2
    $JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$manhattanvcf 
    $JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$euclideanvcf 
    
    echo 'Manhattan VCF: '$manhattanvcf
    echo 'Euclidean VCF: '$euclideanvcf
    
    manhattantsv=$BINDIR/'jasmine_manhattan'.merged.tsv
    euclideantsv=$BINDIR/'jasmine_euclidean'.merged.tsv
    javac $SRCPATH/*.java
    java -cp $SRCPATH VcfToTsv $manhattanvcf $manhattantsv
    echo 'Manhattan TSV: '$manhattantsv
    
    java -cp $SRCPATH VcfToTsv $euclideanvcf $euclideantsv
    echo 'Euclidean TSV: '$euclideantsv
done

