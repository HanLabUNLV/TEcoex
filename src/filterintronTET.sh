
file=$1
filename=$(basename $file)
dir=$(dirname $file)

grep -v 'intron' $file > "${filename}.tmp"
grep ':' "${filename}.tmp" > "${dir}/${filename}2"
rm "${filename}.tmp"
