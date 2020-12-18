BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
ALIGNER='ngmlr'
JASMINESRCPATH='/home/mkirsche/eclipse-workspace/Jasmine/src'
JASMINEEVALSRCPATH='/home/mkirsche/git/JasmineResults/src'
IRISSRCPATH='home/mkirsche/git/Iris/src'

mds=( 5 10 15 20 25 30 35 40 50 60 70 80 90 100 250 500 750 1000 2000 )

for md in ${mds[@]}
do
    echo $md
    HG002FILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/hifi/'$ALIGNER'/'*'_'$md'md_'*`
    HG003FILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg003/hifi/'$ALIGNER'/'*'_'$md'md_'*`
    HG004FILE=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg004/hifi/'$ALIGNER'/'*'_'$md'md_'*`
    echo 'HG002 variants: '$HG002FILE
    echo 'HG003 variants: '$HG003FILE
    echo 'HG004 variants: '$HG004FILE
    
    prefix='jasmine_md'$md
    filelist=$BINDIR/$prefix.filelist.txt
    echo $HG002FILE > $filelist
    echo $HG003FILE >> $filelist
    echo $HG004FILE >> $filelist

    echo 'File list: '$filelist
    
    mergedvcf=$BINDIR/$prefix.merged.vcf

    /home/mkirsche/eclipse-workspace/Jasmine/jasmine file_list=$filelist out_file=$mergedvcf 
    /home/mkirsche/eclipse-workspace/Jasmine/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
    
    echo 'Merged VCF: '$mergedvcf
    
    mergedtsv=$BINDIR/$prefix.merged.tsv
    java -cp $JASMINEEVALSRCPATH VcfToTsv $mergedvcf $mergedtsv
    echo 'Merged TSV: '$mergedtsv
done
