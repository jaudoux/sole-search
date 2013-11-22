more $1 | grep "chr" | awk '{OFS="\t"; print $11, $12, $13}' > temp_$1
more temp_$1 | grep "F" | awk '{OFS="\t"; print $1, $2, $2+200}' >plus_$1
more temp_$1 | grep "R" | awk '{OFS="\t"; print $1, $2-168, $2+32}' >minus_$1
cat plus_$1 minus_$1 | sort | uniq | sort -k 1,1 -k 2n > $2\_full-length.txt
auk '{OFS="\t"; print $1, $2, $2+32}' temp_$1 >$2\_tags.txt
rm temp_$1
rm plus_$1
rm minus_$1






