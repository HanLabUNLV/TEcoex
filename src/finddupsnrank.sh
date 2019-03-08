resultdir=$1

 awk '{print $2}' $resultdir/coexpressed_plus.txt | sort > $resultdir/sorted.plus
 awk '{print $2}' $resultdir/coexpressed_plus.txt | sort | uniq > $resultdir/uniq.plus
 diff $resultdir/sorted.plus $resultdir/uniq.plus | grep '<' | uniq  | sed 's/< //' > $resultdir/duplicate.plus.tmp
 awk -F"|" '{print $1"\t"$2}' $resultdir/duplicate.plus.tmp > $resultdir/duplicate.plus.txt
 awk -F"|" '{print $2}' $resultdir/duplicate.plus.tmp > $resultdir/duplicate.plus.entrez.txt

 awk '{print $2}' $resultdir/coexpressed_minus.txt | sort > $resultdir/sorted.minus
 awk '{print $2}' $resultdir/coexpressed_minus.txt | sort | uniq > $resultdir/uniq.minus
 diff $resultdir/sorted.minus $resultdir/uniq.minus | grep '<' | uniq  | sed 's/< //'  > $resultdir/duplicate.minus.tmp
 awk -F"|" '{print $1"\t"$2}' $resultdir/duplicate.minus.tmp > $resultdir/duplicate.minus.txt
 awk -F"|" '{print $2}' $resultdir/duplicate.minus.tmp > $resultdir/duplicate.minus.entrez.txt

rm -f $resultdir/*.plus $resultdir/*.minus
sort $resultdir/duplicate.plus.tmp > /tmp/dup.plus
sort -k 2b,2 $resultdir/coexpressed_plus.txt > /tmp/coex.plus 
join -1 1 -2 2 /tmp/dup.plus /tmp/coex.plus > $resultdir/coexpressed_plus.duplicates.txt
sort $resultdir/duplicate.minus.tmp > /tmp/dup.minus
sort -k 2b,2 $resultdir/coexpressed_minus.txt > /tmp/coex.minus 
join -1 1 -2 2 /tmp/dup.minus /tmp/coex.minus > $resultdir/coexpressed_minus.duplicates.txt


awk '{print $2"\t"$8}' $resultdir/coexpressed_plus.txt | sort -k 2 -r | awk -F"	" '!_[$1]++' > $resultdir/coexpressed_plus.nr.txt
awk '{print $2"\t-"$8}' $resultdir/coexpressed_minus.txt | sort -k 2 -n | awk -F"	" '!_[$1]++' | sort -n -k 2 -r > $resultdir/coexpressed_minus.nr.txt
cat $resultdir/coexpressed_plus.nr.txt $resultdir/coexpressed_minus.nr.txt | grep -v qval | sort -k 2 -n | awk -F"	" '!_[$1]++' > $resultdir/coexpressed.rankedbyrsquared.rnk
#awk -F'[|\t]' '{print $1"\t"$3}' coexpressed.rankedbyrsquared.tmp > $resultdir/coexpressed.rankedbyrsquared.rnk

sort -k 1b,1 -t"	" $resultdir/coexpressed.rankedbyrsquared.rnk | join -t"	" $resultdir/duplicate.minus.tmp - > $resultdir/duplicate.minus.rsquared.tmp
sort -k 1b,1 -t"	" $resultdir/coexpressed.rankedbyrsquared.rnk | join -t"	" $resultdir/duplicate.plus.tmp - > $resultdir/duplicate.plus.rsquared.tmp
cat $resultdir/duplicate.minus.rsquared.tmp $resultdir/duplicate.plus.rsquared.tmp | sort -k 2 -n -t"	" | awk -F"	" '!_[$1]++' > $resultdir/duplicate.rankedbyrsquared.rnk.tmp
awk -F"|" '{print $1"\t"$2}' $resultdir/duplicate.minus.rsquared.tmp | grep -v ':'  > $resultdir/duplicate.minus.rsquared
awk -F"|" '{print $1"\t"$2}' $resultdir/duplicate.plus.rsquared.tmp | grep -v ':'  > $resultdir/duplicate.plus.rsquared
awk -F"|" '{print $1"\t"$2}' $resultdir/duplicate.rankedbyrsquared.rnk.tmp | grep -v ':'  > $resultdir/duplicate.rankedbyrsquared.rnk
#rm -f *.tmp
