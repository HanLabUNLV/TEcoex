
combine_info() {
inputfile=$1
resultsdir=$2
filename=$(basename "$inputfile")
extension="${filename##*.}"
filename="${filename%.*}"
mkdir "$resultdir/FullInfo.noovp.1k/$filename"

while read LINE
do 
  file="$resultdir/VSTcnts.noovp.1k/$LINE.txt"
  if [ -f "$file" ]; then
    out="$resultdir/FullInfo.noovp.1k/$filename/$LINE.txt"
    awk '{print $2}' $file  > /tmp/tmpfile
    paste ../../data/patient.info /tmp/tmpfile > $out 
  fi
done < $inputfile
}


resultdir='result.rsem.TET'
if [ $# -eq 1 ]; then
  resultdir=$1
fi
mkdir "$resultdir/FullInfo.noovp.1k"
combine_info ../../data/ZNF/DNAlist.txt $resultdir
combine_info ../../data/ZNF/LINElist.txt $resultdir
combine_info ../../data/ZNF/LTRlist.txt $resultdir
combine_info ../../data/ZNF/SINElist.txt $resultdir
combine_info ../../data/ZNF/SVAlist.txt $resultdir

find $resultdir/FullInfo.noovp.1k/ -name '*.txt' | sed 's/^\.\.\/\.\.\///' > $resultdir/FullInfo.noovp.1k/FullInfo.noovp.1k.TElist.txt
