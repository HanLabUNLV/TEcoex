buffer=0
if [ $# ] 
  then 
    buffer=$1
fi


awk -F";" '{print $2"\t"$0}' ~/storage/TETranscript/gtf/hg19/rmsk_TE.sorted.gtf | sed 's/ transcript_id //' | cut -d "	" -f 1-9 | sed 's/"//g'  > /tmp/rmsk.txt
sort -k 1,1 -t"	" /tmp/rmsk.txt | uniq > /tmp/rmsk.sort
awk -F" " '{print $1}' ../../result.rsem.TET.instance/ZNF/significant_instances.txt | awk -F":" '{print $1}' > /tmp/sig.TE.txt
sort -k 1,1 -t"	" /tmp/sig.TE.txt | uniq > /tmp/sig.TE.sort
join -t"	" /tmp/sig.TE.sort /tmp/rmsk.sort > /tmp/sig.TE.pos
awk -F"	" '{print $2"\t"$5"\t"$6"\t"$1"\t"$7"\t"$8}' /tmp/sig.TE.pos | sort -t"	" -k1,1 -k2,2n > /tmp/sig.TE.bed

awk -F";" '{print $1}' ~/storage/TETranscript/gtf/hg19/hg19refGene.gtf  | awk -F"	" '{print $9"\t"$0}' |  cut -d "	" -f 1-9 | grep transcript | sed 's/"//g' | sed 's/gene_id //' > /tmp/refGene.txt
sort -k 1,1 -t"	" /tmp/refGene.txt | uniq > /tmp/refGene.sort
awk -F" " '{print $2}' ../../result.rsem.TET.instance/ZNF/significant_instances.txt | awk -F"|" '{print $1}' > /tmp/sig.gene.txt
sort -k 1,1 -t"	" /tmp/sig.gene.txt | uniq > /tmp/sig.gene.sort
join -t"	" /tmp/sig.gene.sort /tmp/refGene.sort > /tmp/sig.gene.pos
awk -F"	" -v buf="${buffer}" '{if ($5+0 < buf+0) {beg=1} else {beg=0+$5-buf;} print $2"\t"beg"\t"0+$6+buf"\t"$1"\t"$7"\t"$8}' /tmp/sig.gene.pos | sort -t"	" -k1,1 -k2,2n > /tmp/sig.gene.bed

bedtools intersect -wa -wb -a /tmp/sig.TE.bed -b /tmp/sig.gene.bed -sorted > ../../result.rsem.TET.instance/ZNF/sig.TE.ZNF.$buffer.ovp
