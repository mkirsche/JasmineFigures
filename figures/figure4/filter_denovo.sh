BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")

vcf=$BINDIR/denovo.candidates.vcf

java -cp $BINDIR/../../src/ VcfToTsv $vcf

tsv=$BINDIR/denovo.candidates.tsv


