BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
MAX_DIST=50
ALIGNER='ngmlr'
JASMINESRCPATH='/home/mkirsche/eclipse-workspace/Jasmine/src'
JASMINEEVALSRCPATH='/home/mkirsche/git/JasmineResults/src'
IRISSRCPATH='home/mkirsche/git/Iris/src'

PBFILE=`ls '/home/mkirsche/jasmine_data/schatz_pipeline/ash_trio/hg002/pb/'$ALIGNER'/'*'_'$MAX_DIST'md_'*`
ONTFILE=`ls '/home/mkirsche/jasmine_data/schatz_pipeline/ash_trio/hg002/ont/'$ALIGNER'/'*'_'$MAX_DIST'md_'*`
HIFIFILE=`ls '/home/mkirsche/jasmine_data/schatz_pipeline/ash_trio/hg002/hifi/'$ALIGNER'/'*'_'$MAX_DIST'md_'*`

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

/home/mkirsche/eclipse-workspace/Jasmine/jasmine file_list=$filelist out_file=$mergedvcf 
/home/mkirsche/eclipse-workspace/Jasmine/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 
echo 'Merged VCF: '$mergedvcf

mergedtsv=$BINDIR/$prefix.merged.tsv
java -cp $JASMINEEVALSRCPATH VcfToTsv $mergedvcf $mergedtsv
echo 'Merged TSV: '$mergedtsv
