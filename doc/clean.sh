#!/bin/bash

to_del=( *.aux \
		 *.fls \
		 *.log \
		 *.blg \
		 *.bbl \
		 *.fdb_latexmk \
		 *.out \
		 *.pages \
		 *.toc \
		 *.end \
		 *.xdv \
		 *synctex.gz \
		 *Notes.bib )

i=0
while [ $i -lt ${#to_del[@]} ]
do
	if [ -e ${to_del[i]} ]
		then
			echo 'deleting '${to_del[i]}
			rm ${to_del[i]}
	fi;
	i=`expr $i + 1`
done

