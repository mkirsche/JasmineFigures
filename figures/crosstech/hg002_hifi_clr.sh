BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
HIFI_FILE=$BINDIR/../../schatz_pipeline/ash_trio/hg002/hifi/winnowmap/HG002vGRCh38_wm_50md_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf
CLR_FILE=$BINDIR/../../schatz_pipeline/ash_trio/hg002/pb/winnowmap/HG002vGRCh38_wm_50md_PB_sniffles.s2l20.refined.nSVtypes.ism.vcf

JASMINE_SRC_PATH=$BINDIR/../../Jasmine/src
IRIS_SRC_PATH=$BINDIR/../../Jasmine/Iris/src
SRC_PATH=$BINDIR/../../src/

dist=$1

javac $IRIS_SRC_PATH/*.java
javac -cp $IRIS_SRC_PATH $JASMINE_SRC_PATH/*.java
javac $SRC_PATH/*.java

FILELIST=$BINDIR/filelist.txt
echo $HIFI_FILE > $FILELIST
echo $CLR_FILE >> $FILELIST

DISTLIST=$BINDIR/distlist_$i'.txt'
echo $i > $DISTLIST
echo $i >> $DISTLIST

OUTVCF=$BINDIR/merged_$dist'.vcf'

java -cp $IRIS_SRC_PATH:$JASMINE_SRC_PATH Main file_list=$FILELIST out_file=$OUTVCF max_dist=$dist sample_dists=$DISTLIST
java -cp $IRIS_SRC_PATH:$JASMINE_SRC_PATH Main out_file=$OUTVCF --postprocess_only --dup_to_ins
java -cp $SRC_PATH VcfToTsv $OUTVCF
