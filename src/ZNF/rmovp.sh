if [ $# ] 
  then 
    SIGFILE=$1
    OVPFILE=$2
fi

#echo $SIGFILE
#echo $OVPFILE

while read LINE
do 
  TE=$(echo $LINE | awk '{print $1}') 
  ZNF=$(echo $LINE | awk '{print $2}') 
  TE=$(echo ${TE} | awk -F":" '{print $1}')
  ZNF=$( echo "${ZNF}"| awk -F"|" '{print $1}')

  found=0
  grepcmd="grep -w ${TE} $OVPFILE "
  check_ovp=`$grepcmd`
  if [[ $check_ovp ]]; then 
    if echo "$check_ovp" | grep -wq "${ZNF}"; then
      found=1
    fi 
  fi
  if [ $found -eq 0 ]; then
    echo $LINE
  fi
done < $SIGFILE 

