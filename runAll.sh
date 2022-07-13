for FILE in ./NewSumTab/*.tsv; 
do echo $FILE; 
cat $FILE;
./main $FILE 0;
break
done