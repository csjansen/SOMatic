#!/bin/bash
#$ -N MM3A
#$ -q sam128,sam,bio
#$ -pe openmp 1

rm "$3".matrix 2> /dev/null
awk '{ a[FNR] = (a[FNR] ? a[FNR] FS : "") $6 } END { for(i=2;i<=FNR;i++) print a[i] }' $(ls -1v $1) > "$3".matrix

#to create a sample.list with the Cell names from the Matrix
rm $2 2> /dev/null
for i in `ls -1v $1`; do echo ${i%.*} >> $2; done;
sed -i -r "s/(.*)\..*\/([^\/\.]*)\..*/\1.\2/g" $2

sed -i 's/ /\t/g' "$3".matrix
topfile=`ls -1v $1 | head -1`
cut -f1 $topfile | tail -n +2 > "$3".genes
paste "$3".genes "$3".matrix > $3
rm "$3".genes 2> /dev/null
rm "$3".matrix 2> /dev/null
