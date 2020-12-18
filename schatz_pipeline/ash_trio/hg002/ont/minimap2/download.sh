BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp mkirsche@jhu.edu@gateway2.marcc.jhu.edu:/work-zfs/mschatz1/giab/mm2_AshkenazimTrio/sniffles_md_parameters/vmds/svs/refined/*HG002_2020*ONT*.refined.nSVtypes.ism.vcf $BINDIR
