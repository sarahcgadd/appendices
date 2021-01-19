#!/bin/bash
CORES=1

echo "FILEREF;BOOTS;NODES;PARAMREF" > /home/home02/mt15scg/patternsim2/fileref_pattern.txt

reference=1

for BOOTS in 1000
do

for NODES in  20
do

for PARAM in {1..8} 
do

for SEED in {1..2}
do


echo "paramref <- $PARAM" > bigsim_pattern_$reference.R
echo "fileref<-$reference; nstations<- $NODES; nsims<-1; boots<-$BOOTS; cores<-$CORES" >> bigsim_pattern_$reference.R
echo "set.seed($SEED)" >> bigsim_pattern_$reference.R
cat bigsim_pattern.R >> bigsim_pattern_$reference.R
echo "write.csv(result, paste(\"/home/home02/mt15scg/patternsim2/results_pattern_$reference.csv\", sep=\"\"))" >> bigsim_pattern_$reference.R

echo "$reference;$BOOTS;$NODES;$PARAM" >> OUTPUTFOLDER/fileref_pattern.txt

chmod 744 bigsim_pattern_$reference.R

reference=$((reference+1))

done

done

done

done

reference=$((reference-1))

echo "#!/bin/bash" > test_array_2.sh

echo "#$ -pe smp $CORES" >> test_array_2.sh

echo "#$ -t 1-$reference" >> test_array_2.sh

cat gen_subscript.sh >> test_array_2.sh

chmod 744 test_array_2.sh

echo "#!/bin/bash" > collate_results_pattern.sh

echo "cp results_pattern_1.csv allresults_pattern.csv" >> collate_results_pattern.sh

echo "for i in {2..$reference}" >> collate_results_pattern.sh

cat collate_pattern.sh >> collate_results_pattern.sh

chmod 744 collate_results_pattern.sh

cp collate_results_pattern.sh OUTPUTFOLDER/collate_results_pattern.sh
