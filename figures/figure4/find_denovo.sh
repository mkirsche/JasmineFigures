if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

MAX_DIST=50
ALIGNER=ngmlr
JASMINEEVALSRCPATH='/home/mkirsche/git/JasmineResults/src'

hifi_child_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/hifi/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
hifi_father_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg003/hifi/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
hifi_mother_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg004/hifi/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
clr_child_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/pb/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
clr_father_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg003/pb/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
clr_mother_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg004/pb/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
ont_child_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg002/ont/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
ont_father_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg003/ont/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`
ont_mother_vcf=`ls $BINDIR'/../../schatz_pipeline/ash_trio/hg004/ont/'$ALIGNER'/'*'_'$MAX_DIST'md_'*ism.vcf`

prefix=denovo

filelist=$BINDIR/$prefix.filelist.txt
echo $hifi_child_vcf > $filelist
echo $hifi_father_vcf >> $filelist
echo $hifi_mother_vcf >> $filelist
echo $clr_child_vcf >> $filelist
echo $clr_father_vcf >> $filelist
echo $clr_mother_vcf >> $filelist
echo $ont_child_vcf >> $filelist
echo $ont_father_vcf >> $filelist
echo $ont_mother_vcf >> $filelist

mergedvcf=$BINDIR/$prefix.merged.vcf

/home/mkirsche/eclipse-workspace/Jasmine/jasmine file_list=$filelist out_file=$mergedvcf 
/home/mkirsche/eclipse-workspace/Jasmine/jasmine --dup_to_ins --postprocess_only out_file=$mergedvcf 

echo 'Merged VCF: '$mergedvcf

mergedtsv=$BINDIR/$prefix.merged.tsv
java -cp $JASMINEEVALSRCPATH VcfToTsv $mergedvcf $mergedtsv

Rscript $BINDIR/../jasmine_plot/find_denovo_candidates.R

candidatesvcf=$BINDIR/$prefix.candidates.vcf
cat $mergedvcf | grep -v 'IMPRECISE;' | grep 'IS_SPECIFIC=1' | grep 'SUPP_VEC=100100100' > $candidatesvcf