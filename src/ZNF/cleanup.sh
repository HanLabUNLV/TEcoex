inputfile=$1

lookupdir="../../result.rsem.TET.instance/ZNF.2rm"
mvdir="../../result.rsem.TET.instance/ZNF"
while IFS=$'\t' read -r column1 column2 column3 columnrest
do 
  filename1="${column3}/${column2}_${column1}"
  filename2="${column2}_${column1}_${column3}"
  #echo "${lookupdir}/${filename1}_results_positive_lm.txt"
  #echo "${lookupdir}/${filename2}_coexpressed_plus.txt"
  #echo "${lookupdir}/${filename2}_coexpressed_minus.txt" 

  if [ -f "${lookupdir}/${filename1}_results_positive_lm.txt" ]; then 
    cmd="mv ${lookupdir}/${filename1}_results_*.txt ${mvdir}/${column3}/"
    `$cmd`
    echo "${column1}	${column2}	${column3}"
  fi
  #if [ ! -f "${lookupdir}/${filename2}_coexpressed_plus.txt" ]; then 
    #cmd="mv ${lookupdir}/${filename2}_coexpressed_plus.txt ${mvdir}"
    #`$cmd`
  #fi
  #if [ ! -f "${lookupdir}/${filename2}_coexpressed_minus.txt" ]; then 
    #cmd="mv ${lookupdir}/${filename2}_coexpressed_minus.txt ${mvdir}"
    #`$cmd`
  #fi
done < $inputfile

