#!/bin/bash

for FILEREF in {1..1000}
do
	for CHAINNUM in {1..4}
	do
	NEWNUM=$(( 4*(FILEREF-1)+CHAINNUM ))
	echo "load(file='renvironment0_$FILEREF.RData')" > singlesim_1_$NEWNUM.R
	echo "fileref<-$FILEREF" >> singlesim_1_$NEWNUM.R
	echo "set.seed($FILEREF)" >> singlesim_1_$NEWNUM.R
	echo "chainnum<-$CHAINNUM" >> singlesim_1_$NEWNUM.R
	cat singlesim_1.R >> singlesim_1_$NEWNUM.R
	#chmod 744 singlesim_1_$NEWNUM.R
	done
	
	echo "fileref<-$FILEREF" > singlesim_2_$FILEREF.R
        cat singlesim_2.R >> singlesim_2_$FILEREF.R
        chmod 744 singlesim_2_$FILEREF.R

	echo "fileref<-$FILEREF" > singlesim_0_$FILEREF.R
        cat singlesim_0.R >> singlesim_0_$FILEREF.R
        chmod 744 singlesim_0_$FILEREF.R

	echo "fileref<-$FILEREF" > singlesim_4_$FILEREF.R
	cat singlesim_4.R >> singlesim_4_$FILEREF.R
	echo "write.csv(result, paste(\"OUTPUTFOLDER/results_pattern_$FILEREF.csv\", sep=\"\"))" >> singlesim_4_$FILEREF.R
	chmod 744 singlesim_4_$FILEREF.R

	for LOOPNUM in {1..40}
	do

		filenum=$(( 40*(FILEREF-1)+LOOPNUM ))
		echo "fileref<-$FILEREF" > singlesim_3_$filenum.R
		echo "load(paste('renvironment_2_',fileref,'.RData',sep=''))" >> singlesim_3_$filenum.R
		echo "loopnum<-$LOOPNUM" >> singlesim_3_$filenum.R
		cat singlesim_3.R >> singlesim_3_$filenum.R
		chmod 744 singlesim_3_$filenum.R

	done

done


echo "#!/bin/bash" > sim_array1.sh

echo "#$ -pe smp 1" >> sim_array1.sh

echo "#$ -t 1-$NEWNUM" >> sim_array1.sh

cat gen_subscript1.sh >> sim_array1.sh

chmod 744 sim_array1.sh



echo "#!/bin/bash" > sim_array3.sh

echo "#$ -pe smp 1" >> sim_array3.sh

echo "#$ -t 1-$filenum" >> sim_array3.sh

cat gen_subscript3.sh >> sim_array3.sh

chmod 744 sim_array3.sh


echo "#!/bin/bash" > sim_array2.sh

echo "#$ -pe smp 1" >> sim_array2.sh

echo "#$ -t 1-$FILEREF" >> sim_array2.sh

cat gen_subscript2.sh >> sim_array2.sh

chmod 744 sim_array2.sh


echo "#!/bin/bash" > sim_array0.sh

echo "#$ -pe smp 1" >> sim_array0.sh

echo "#$ -t 1-$FILEREF" >> sim_array0.sh

cat gen_subscript0.sh >> sim_array0.sh

chmod 744 sim_array0.sh


echo "#!/bin/bash" > sim_array4.sh

echo "#$ -pe smp 1" >> sim_array4.sh

echo "#$ -t 1-$FILEREF" >> sim_array4.sh

cat gen_subscript4.sh >> sim_array4.sh

chmod 744 sim_array4.sh


echo "#!/bin/bash" > collate_results_pattern.sh

echo "cp results_pattern_1.csv allresults_pattern.csv" >> collate_results_pattern.sh

echo "for i in {2..$FILEREF}" >> collate_results_pattern.sh

cat collate_pattern.sh >> collate_results_pattern.sh

chmod 744 collate_results_pattern.sh

cp collate_results_pattern.sh /OUTPUTFOLDER/collate_results_pattern.sh

   
