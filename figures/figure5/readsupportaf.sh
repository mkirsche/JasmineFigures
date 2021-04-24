BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")

JASMINEPATH=$BINDIR'/../../Jasmine'
echo 'Compiling'
javac -cp $JASMINEPATH/Iris/src:$JASMINEPATH/src /home/mkirsche/eclipse-workspace/SvPopulationAnalysis/src/*.java
echo 'Adding Read Support'
java -cp $JASMINEPATH/Iris/src:$JASMINEPATH/src:/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/src AddReadSupport table_file=jasmine_pop_final_specprec.txt out_file=all4_pop_jasmine_support.txt vcf_filelist=$BINDIR/rockfishfilelist.txt
java -cp $JASMINEPATH/Iris/src:$JASMINEPATH/src:/home/mkirsche/eclipse-workspace/SvPopulationAnalysis/src FilterSpecific
