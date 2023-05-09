#!/bin/bash

#SBATCH --mem=10G

filename="uniid.txt"
while read -r line; do
	name="$line"
	wget https://www.uniprot.org/uniprot/"$name".txt
	cat "$name".txt | grep "^DR\s*GO;\s*GO:[0-9]*;\s*P" > try.txt
	sed -i -e "s/^/$name;/" try.txt
	cat try.txt >> out.txt
	rm "$name".txt
done < "$filename"
