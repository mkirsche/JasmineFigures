BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")

git submodule update --init --recursive

cd $BINDIR/SURVIVOR/Debug
make

cd $BINDIR/Jasmine
./build_jar.sh
javac -cp Iris/src src/*.java

cd $BINDIR/src
javac *.java

cd $BINDIR/figures/figure2/SvPopulationAnalysis
javac -cp $BINDIR/Jasmine/src src/*.java
