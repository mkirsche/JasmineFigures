BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
for i in `ls $BINDIR/*/download.sh`
do
  echo $i
  $i
done
