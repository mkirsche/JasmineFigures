BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp mkirsche@jhu.edu@gateway2.marcc.jhu.edu:/scratch/groups/mschatz1/aganezov/CHM13/variants/noAD_GRCh38_ref/GIABvGRCh38_wm/svs/refined/HG004vGRCh38_wm_PBCCS_sniffles.s2l20.refined.nSVtypes.ism.vcf $BINDIR
