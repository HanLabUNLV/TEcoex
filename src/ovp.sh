buffer=0
if [ $# ] 
  then 
    TEdatadir=$1
    buffer=$2
fi



awk -F";" '{print $2"\t"$0}' ~/storage/TETranscript/gtf/hg19/rmsk_TE.sorted.gtf | sed 's/ transcript_id //' | cut -d "	" -f 1-9 | sed 's/"//g'  > /tmp/rmsk.txt
sort -k 1,1 -t"	" /tmp/rmsk.txt | uniq > /tmp/rmsk.sort
awk -F" " '{print $1}' ../data/${TEdatadir}/rawcnts.new.maxgt5.txt | awk -F":" '{print $1}' > /tmp/exp.TE.txt
sort -k 1,1 -t"	" /tmp/exp.TE.txt | uniq > /tmp/exp.TE.sort
join -t"	" /tmp/exp.TE.sort /tmp/rmsk.sort > /tmp/exp.TE.pos
awk -F"	" '{print $2"\t"$5"\t"$6"\t"$1"\t"$7"\t"$8}' /tmp/exp.TE.pos | sort -t"	" -k1,1 -k2,2n > /tmp/exp.TE.bed

awk -F";" '{print $1}' ~/storage/TETranscript/gtf/hg19/hg19refGene.gtf  | awk -F"	" '{print $9"\t"$0}' |  cut -d "	" -f 1-9 | grep transcript | sed 's/"//g' | sed 's/gene_id //' > /tmp/refGene.txt
awk -F"	" -v buf="${buffer}" '{if ($5+0 < buf+0) {beg=1} else {beg=0+$5-buf;} print $2"\t"beg"\t"0+$6+buf"\t"$1"\t"$7"\t"$8}' /tmp/refGene.txt | sort -t"	" -k1,1 -k2,2n > /tmp/all.gene.bed


bedtools intersect -wa -v -a /tmp/exp.TE.bed -b /tmp/all.gene.bed -sorted > ../result.rsem.${TEdatadir}.instance/exp.TE.gene.$buffer.noovp
awk '{print $1"\t"$2"\t"$3"\t"$4}' ../result.rsem.${TEdatadir}.instance/exp.TE.gene.$buffer.noovp | sort | uniq > ../result.rsem.${TEdatadir}.instance/gene.noovp.TEs.bed
bedtools intersect -wa -a ../result.rsem.${TEdatadir}.instance/gene.noovp.TEs.bed -b ~/storage/TETranscript/gtf/hg19/hg19_rmsk_uniq.bed | sort | uniq > ../result.rsem.${TEdatadir}.instance/gene.noovp.uniq.TEs.bed


