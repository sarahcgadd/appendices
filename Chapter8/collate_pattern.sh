do
cat results_pattern_$i.csv | sed '1d' >> allresults_pattern.csv

done
