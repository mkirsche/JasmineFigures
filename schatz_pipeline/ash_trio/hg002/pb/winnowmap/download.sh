BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp mkirsch5@login.rockfish.jhu.edu:/scratch4/mschatz1/Jasmine/analysis/ASHKENAZIvGRCh38_wm/svs/refined/HG002vGRCh38_wm_PB_sniffles.s2l20.refined.nSVtypes.ism.vcf $BINDIR
