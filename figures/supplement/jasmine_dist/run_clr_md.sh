BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
ALIGNER='ngmlr'
JASMINEPATH=$BINDIR'/../../../Jasmine'
SRCPATH=$BINDIR'/../../../src'

mds=( 5 10 20 30 40 50 60 70 80 90 100 250 500 1000 )

for md in ${mds[@]}
do
    echo $md
    HG002FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg002/pb/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    HG003FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg003/pb/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    HG004FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg004/pb/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    echo 'HG002 variants: '$HG002FILE
    echo 'HG003 variants: '$HG003FILE
    echo 'HG004 variants: '$HG004FILE
    
    prefix='jasmine_md_clr'$md
    filelist=$BINDIR/$prefix.filelist.txt
    echo $HG002FILE > $filelist
    echo $HG003FILE >> $filelist
    echo $HG004FILE >> $filelist

    echo 'File list: '$filelist
    
    mergedvcf=$BINDIR/$prefix.merged.vcf

    #$JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf 
    #$JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
    
    echo 'Merged VCF: '$mergedvcf
    
    mergedtsv=$BINDIR/$prefix.merged.tsv
    javac $SRCPATH/*.java
    java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv
    echo 'Merged TSV: '$mergedtsv
done

for md in ${mds[@]}
do
    echo $md
    HG002FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg002/ont/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    HG003FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg003/ont/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    HG004FILE=`ls $BINDIR'/../../../schatz_pipeline/ash_trio/hg004/ont/'$ALIGNER'/'*'_'$md'md_'*.ism.vcf`
    echo 'HG002 variants: '$HG002FILE
    echo 'HG003 variants: '$HG003FILE
    echo 'HG004 variants: '$HG004FILE
    
    prefix='jasmine_md_ont'$md
    filelist=$BINDIR/$prefix.filelist.txt
    echo $HG002FILE > $filelist
    echo $HG003FILE >> $filelist
    echo $HG004FILE >> $filelist

    echo 'File list: '$filelist
    
    mergedvcf=$BINDIR/$prefix.merged.vcf

    #$JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf --output_genotypes
    #$JASMINEPATH/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
    
    echo 'Merged VCF: '$mergedvcf
    
    mergedtsv=$BINDIR/$prefix.merged.tsv
    javac $SRCPATH/*.java
    java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv
    echo 'Merged TSV: '$mergedtsv
done

Rscript $BINDIR/../../jasmine_plot/jasmine_dist_discordance.R

