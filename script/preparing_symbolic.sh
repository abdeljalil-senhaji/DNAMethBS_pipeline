#!/usr/bin/env bash

# $1 = Dossier des échantillons
# $2 = Fichier des résultats


for file in $(ls $1)
do
	# if is pattern is good just create simple link
	echo $file
	if echo $file | grep -q '_R1.fastq.gz\|_R2.fastq.gz' ; then   		
		ln -s $1/$file $2/$file 
	elif echo $file | grep -Eq '(_[0-9]*_)1_(.*)\.(fq|FQ|FASTQ|fastq)\.gz' ; then 
		echo toto
		newName=$(echo $file | sed -r "s/(_[0-9]*_)1_(.*)\.(fq|FQ|FASTQ|fastq)\.gz/\1\2_R1.fastq.gz/g")  		
		ln -s $1/$file $2/$newName
	elif echo $file | grep -Eq '(_[0-9]*_)2_(.*)\.(fq|FQ|FASTQ|fastq)\.gz' ; then 
		echo tata
		newName=$(echo $file | sed -r "s/(_[0-9]*_)2_(.*)\.(fq|FQ|FASTQ|fastq)\.gz/\1\2_R2.fastq.gz/g")  		
		ln -s $1/$file $2/$newName
	elif echo $file | grep -Eq '(.*)(_1|^1[_-]|_[Rr]1|^[Rr]1[_-])(.*)\.(fq|FQ|FASTQ|fastq)\.gz' ; then 
		echo tata
		newName=$(echo $file | sed -r "s/(.*)(_1|^1[_-]|_[Rr]1|^[Rr]1[_-])(.*)\.(fq|FQ|FASTQ|fastq)\.gz/\1\3_R1.fastq.gz/g")  		
		ln -s $1/$file $2/$newName
	elif echo $file | grep -Eq '(.*)(_2|^2[_-]|_[Rr]2|^[Rr]2[_-])(.*)\.(fq|FQ|FASTQ|fastq)\.gz' ; then 
		echo tata
		newName=$(echo $file | sed -r "s/(.*)(_2|^2[_-]|_[Rr]2|^[Rr]2[_-])(.*)\.(fq|FQ|FASTQ|fastq)\.gz/\1\3_R2.fastq.gz/g")  		
		ln -s $1/$file $2/$newName
	fi
done



