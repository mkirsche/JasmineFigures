BINDIR=$(dirname "$(readlink -f "$0" || echo "$(echo "$0" | sed -e 's,\\,/,g')")")
scp mkirsch5@login.rockfish.jhu.edu:/scratch4/mschatz1/mkirsche/JasmineFigures/figures/figure5/all4_pop.*augmented.txt $BINDIR
cat $BINDIR/all4_pop.jasmine_augmented.txt | awk '{if(($44 == "1" || $44 == "SPECIFIC_FLAG") && ($45 == "1" || $45 == "PRECISE_FLAG")) {print}}' > $BINDIR/all4_pop.jasmine.specprec.txt
cat $BINDIR/all4_pop.survivor_augmented.txt | awk '{if(($44 == "1" || $44 == "SPECIFIC_FLAG") && ($45 == "1" || $45 == "PRECISE_FLAG")) {print}}' > $BINDIR/all4_pop.survivor.specprec.txt
cat $BINDIR/all4_pop.svtools_augmented.txt | awk '{if(($44 == "1" || $44 == "SPECIFIC_FLAG") && ($45 == "1" || $45 == "PRECISE_FLAG")) {print}}' > $BINDIR/all4_pop.svtools.specprec.txt
cat $BINDIR/all4_pop.svimmer_augmented.txt | awk '{if(($44 == "1" || $44 == "SPECIFIC_FLAG") && ($45 == "1" || $45 == "PRECISE_FLAG")) {print}}' > $BINDIR/all4_pop.svimmer.specprec.txt
cat $BINDIR/all4_pop.svpop_augmented.txt | awk '{if(($44 == "1" || $44 == "SPECIFIC_FLAG") && ($45 == "1" || $45 == "PRECISE_FLAG")) {print}}' > $BINDIR/all4_pop.svpop.specprec.txt

