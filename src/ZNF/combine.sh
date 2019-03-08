
combine_info() {
inputfile=$1
resultsdir=$2
filename=$(basename "$inputfile")
extension="${filename##*.}"
filename="${filename%.*}"
mkdir "$resultdir/FullInfo/$filename"

while read LINE
do 
  file="$resultdir/VSTcnts/$LINE.txt"
  out="$resultdir/FullInfo/$filename/$LINE.txt"
  awk '{print $2}' $file  > /tmp/tmpfile
  paste ../../data/patient.info /tmp/tmpfile > $out 
done < $inputfile
}


resultdir='result.rsem.TET'
if [ $# -eq 1 ]; then
  resultdir=$1
fi
mkdir "$resultdir/FullInfo"
combine_info ../../data/ZNF/DNAlist.txt $resultdir
combine_info ../../data/ZNF/LINElist.txt $resultdir
combine_info ../../data/ZNF/LTRlist.txt $resultdir
combine_info ../../data/ZNF/SINElist.txt $resultdir
combine_info ../../data/ZNF/SVAlist.txt $resultdir

find $resultdir/FullInfo/ -name '*.txt' | sed 's/^\.\.\/\.\.\///' > $resultdir/FullInfo/FullInfo.TElist.txt
