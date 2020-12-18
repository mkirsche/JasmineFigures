BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp mkirsche@jhu.edu@gateway2.marcc.jhu.edu:/work-zfs/mschatz1/giab/AshkenazimTrio/sniffles_params_workdir/vmds/svs/refined/HG004_*PB_*.refined.nSVtypes.ism.vcf $BINDIR
