
datadir=$1
cd $datadir 

awk '{print $1}' ../TET/BRCA/TCGA-A7-A0CE.discount.instance.cntTable > rownames.txt
mkdir cnts


for file in */*discount.instance.cntTable; 
do 
  tmp=$file
  tmp=${tmp##*/}
  cntfilename=${tmp%.discount.instance.cntTable}   # remove suffix starting with "_"
  echo $cntfilename

  awk '{print $2}' $file | tail -4496028 > tmpcnt
  echo $cntfilename > cnts/$cntfilename.txt 
  cat tmpcnt >> cnts/$cntfilename.txt 
done


paste rownames.txt cnts/* > rawcnts.new.txt

awk 'BEGIN{max=0}{for(i=2;i<=NF;i++){max=(max>$i)?max:$i}if(max>5)print $0;max=0}' rawcnts.new.txt > rawcnts.new.maxgt5.txt &
