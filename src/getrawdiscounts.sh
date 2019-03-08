
ele_name=$1
echo ${ele_name}

grep ${ele_name} ../data/TET/*/*ele.cntTable > ${ele_name}.raw.txt
grep -v 'discount' ${ele_name}.raw.txt |  grep -v intron > ../result.rsem.TET/${ele_name}.rawcnts.txt 
grep 'discount' ${ele_name}.raw.txt > ../result.rsem.TET/${ele_name}.discount.txt
rm ${ele_name}.raw.txt 
 
grep ${ele_name} ../data/TET.uniq/*/*ele.cntTable > ${ele_name}.raw.txt
grep -v 'discount' ${ele_name}.raw.txt |  grep -v intron > ../result.rsem.TET.uniq/${ele_name}.rawcnts.txt 
grep 'discount' ${ele_name}.raw.txt > ../result.rsem.TET.uniq/${ele_name}.discount.txt
rm ${ele_name}.raw.txt 



grep ${ele_name} ../data/TET.bowtie/*/*ele.cntTable > ${ele_name}.raw.txt
grep -v 'discount' ${ele_name}.raw.txt |  grep -v intron > ../result.rsem.TET.bowtie/${ele_name}.rawcnts.txt 
grep 'discount' ${ele_name}.raw.txt > ../result.rsem.TET.bowtie/${ele_name}.discount.txt
rm ${ele_name}.raw.txt 
