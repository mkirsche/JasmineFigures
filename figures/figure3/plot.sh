BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
Rscript $BINDIR/../jasmine_plot/crosstech.R
