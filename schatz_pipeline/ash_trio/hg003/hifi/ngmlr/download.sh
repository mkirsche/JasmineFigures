BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp temp-mkirsche@rflogin01.rockfish.jhu.edu:/scratch4/temp-saganez1/CHM13/GRCh381KGP/GIABvGRCh38_ngmlr/GIABvGRCh38_ngmlr_sniffles_mds/svs/refined/HG003*.refined.nSVtypes.ism.vcf $BINDIR
