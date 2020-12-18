if [ "$(uname -s)" = 'Linux' ]; then
    BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
else
    BINDIR=$(dirname "$(readlink "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
fi

SRCPATH=$BINDIR'/../../src'
TABLESRCPATH=$BINDIR'/../figure2/SvPopulationAnalysis/src'
SURVIVOR_PATH=$BINDIR'/../../SURVIVOR'
JASMINEPATH=$BINDIR'/../../Jasmine'
prefix='population'

filelist=$BINDIR/$prefix.filelist.txt
mergedvcf=$BINDIR/$prefix.merged.vcf
$JASMINEPATH/jasmine file_list=$filelist out_file=$mergedvcf

mergedtsv=$BINDIR/$prefix.merged.tsv
java -cp $SRCPATH VcfToTsv $mergedvcf $mergedtsv

jasminesmalltable=$BINDIR/$prefix.jasminesmalltable.txt
jasminetable=$BINDIR/$prefix.jasminetable.txt
java -cp $JASMINEPATH/src:$TABLESRCPATH BuildMergingTable vcf_file=$mergedvcf out_file=$jasminesmalltable vcf_filelist=$filelist
java -cp $JASMINEPATH/src:$TABLESRCPATH AugmentMergingTable table_file=$jasminesmalltable out_file=$jasminetable vcf_filelist=$filelist

survivorvcf=$BINDIR/$prefix.survivor.vcf
$SURVIVOR_PATH/Debug/SURVIVOR merge $filelist 1000 1 1 1 0 1 $survivorvcf
survivorsmalltable=$BINDIR/$prefix.survivorsmalltable.txt
survivortable=$BINDIR/$prefix.survivortable.txt
java -cp $JASMINEPATH/src:$TABLESRCPATH BuildMergingTable vcf_file=$survivorvcf out_file=$survivorsmalltable vcf_filelist=$filelist mode=survivor
java -cp $JASMINEPATH/src:$TABLESRCPATH AugmentMergingTable table_file=$survivorsmalltable out_file=$survivortable vcf_filelist=$filelist


