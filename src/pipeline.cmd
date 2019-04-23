
##################### rsem genes ########################
Rscript L1HSstatistics.r &> rsem.log 
Rscript makeHeatmap.r result.rsem.TET &> heatmap.log
Rscript gene2L1HS.r  &> lm.log 
Rscript gene2oldLINE.r  &> gene2oldLINE.log 

bash finddupsnrank.sh ../result.rsem.TET/gene2L1HS/
Rscript clusterprofiler.r result.rsem.TET/gene2L1HS/
Rscript david.r result.rsem.TET/gene2L1HS/
Rscript RECscore.r ../result.rsem.TET/gene2L1HS/

bash finddupsnrank.sh ../result.rsem.TET/gene2oldLINE/
Rscript clusterprofiler.r result.rsem.TET/gene2oldLINE/
Rscript david.r result.rsem.TET/gene2oldLINE/
Rscript RECscore.r ../result.rsem.TET/gene2oldLINE/

Rscript WGCNA.L1HS.tissueconsensus.r 

Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ BLCA bladder 7 &> wgcna.BLCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ BRCA breast 20 &> wgcna.BRCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ COADREAD colonrectum 10 &> wgcna.COADREAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ ESCASTAD esophagusstomach 12 &> wgcna.ESCASTAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ HNSC headneck 18 &> wgcna.HNSC.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ KICHKIRCKIRP kidney 6 &> wgcna.KICHKIRCKIRP.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ LIHC liver 8 &> wgcna.LICH.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ LUADLUSC lung 12 &> wgcna.LUAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ PRAD prostate 9 &> wgcna.PRAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ THCA thyroid 20 &> wgcna.THCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET/ UCEC endometrium 12 &> wgcna.UCEC.log &

##################### TEtranscript instance ########################
filterinstancecnts.sh ../data/TET
ovp.sh TET 10000
Rscript L1HSstatistics.instance.r &> instance.log 
Rscript makeHeatmap.r result.rsem.TET.instance &> heatmap.instance.log 
Rscript makeHeatmap.r result.rsem.TET.instance noovp &> heatmap.instance.noovp.log 

# run ZNF scripts before splitVSTcnts.r
Rscript splitVSTcnts.r
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.BLCA.sig.txt bladder 7 &> wgcna.instance.BLCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.BRCA.sig.txt breast 14 &> wgcna.instance.BRCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.COADREAD.sig.txt colonrectum 10 &> wgcna.instance.COADREAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.ESCASTAD.sig.txt esophagusstomach 16 &> wgcna.instance.ESCASTAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.HNSC.sig.txt headneck 18 &> wgcna.instance.HNSC.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.KICHKIRCKIRP.sig.txt kidney 6 &> wgcna.instance.KICHKIRCKIRP.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.LIHC.sig.txt liver 9 &> wgcna.instance.LICH.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.LUADLUSC.sig.txt lung 6 &> wgcna.instance.LUAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.PRAD.sig.txt prostate 10 &> wgcna.instance.PRAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.THCA.sig.txt thyroid 7 &> wgcna.instance.THCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.instance/ VSTcnt.UCEC.sig.txt endometrium 16 &> wgcna.instance.UCEC.log &

##################### TEtranscript genes ########################
Rscript L1HSstatistics.r data/TET data/TET &> TET.log 
Rscript gene2L1HS.r result.TET.TET  &> lmTET.log 
Rscript gene2oldLINE.r result.TET.TET  &> lmTET.oldLINE.log 

bash finddupsnrank.sh ../result.TET.TET/gene2L1HS/
Rscript clusterprofiler.r result.TET.TET/gene2L1HS/
Rscript david.r result.TET.TET/gene2L1HS/

bash finddupsnrank.sh ../result.TET.TET/gene2oldLINE/
Rscript clusterprofiler.r result.TET.TET/gene2oldLINE/
Rscript david.r result.TET.TET/gene2oldLINE/


##################### TEtranscript unique ########################
Rscript L1HSstatistics.r data/rsem data/TET.uniq &> TET.uniq.log 
Rscript gene2L1HS.r result.rsem.TET.uniq &> lmuniq.log 
Rscript gene2oldLINE.r result.rsem.TET.uniq  &> lmuniq.oldLINE.log 
Rscript gene2LINE.coef.r result.rsem.TET.uniq/gene2L1HS 0.001
Rscript gene2LINE.coef.r result.rsem.TET.uniq/gene2oldLINE 0.001

bash finddupsnrank.sh ../result.rsem.TET.uniq/gene2L1HS/
Rscript clusterprofiler.r result.rsem.TET.uniq/gene2L1HS/
Rscript david.r result.rsem.TET.uniq/gene2L1HS/

bash finddupsnrank.sh ../result.rsem.TET.uniq/gene2oldLINE/
Rscript clusterprofiler.r result.rsem.TET.uniq/gene2oldLINE/
Rscript david.r result.rsem.TET.uniq/gene2oldLINE/


##################### bowtie ########################
Rscript L1HSstatistics.r data/rsem data/TET.bowtie &> TET.bowtie.log 
Rscript makeHeatmap.r result.rsem.TET.bowtie &> heatmap.bowtie.log 
Rscript gene2L1HS.r result.rsem.TET.bowtie &> lmbowtie.log 
Rscript gene2LINE.coef.r result.rsem.TET.bowtie/gene2L1HS 0.05
Rscript gene2LINE.coef.r result.rsem.TET.bowtie/gene2oldLINE 0.05

bash finddupsnrank.sh ../result.rsem.TET.bowtie/gene2L1HS/
Rscript clusterprofiler.r result.rsem.TET.bowtie/gene2L1HS/
Rscript david.r result.rsem.TET.bowtie/gene2L1HS/

bash finddupsnrank.sh ../result.rsem.TET.bowtie/gene2oldLINE/
Rscript clusterprofiler.r result.rsem.TET.bowtie/gene2oldLINE/
Rscript david.r result.rsem.TET.bowtie/gene2oldLINE/

Rscript WGCNA.L1HS.tissueconsensus.r result.rsem.TET.bowtie
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ BLCA bladder 8 &> wgcna.BLCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ BRCA breast 18 &> wgcna.BRCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ COADREAD colonrectum 10 &> wgcna.COADREAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ ESCASTAD esophagusstomach 12 &> wgcna.ESCASTAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ HNSC headneck 16 &> wgcna.HNSC.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ KICHKIRCKIRP kidney 7 &> wgcna.KICHKIRCKIRP.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ LIHC liver 8 &> wgcna.LICH.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ LUADLUSC lung 12 &> wgcna.LUAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ PRAD prostate 9 &> wgcna.PRAD.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ THCA thyroid 20 &> wgcna.THCA.log &
Rscript WGCNA.L1HS.bytissue.r result.rsem.TET.bowtie/ UCEC endometrium 12 &> wgcna.UCEC.log &

filterinstancecnts.sh ../data/TET.bowtie/
ovp.sh TET.bowtie 10000
Rscript L1HSstatistics.instance.r data/rsem/ data/TET.bowtie/ &> bowtie.instance.log 
Rscript makeHeatmap.r result.rsem.TET.bowtie.instance &> heatmap.bowtie.instance.log 
Rscript makeHeatmap.r result.rsem.TET.bowtie.instance noovp &> heatmap.bowtie.instance.noovp.log 

# run ZNF scripts before splitVSTcnts.r
Rscript splitVSTcnts.r ../result.rsem.TET.bowtie.instance/
#Rscript splitVSTcnts.r ../result.rsem.TET.bowtie.instance/ noovp

Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.BLCA.sig.txt bladder 0 &> wgcna.bowtie.instance.BLCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.BRCA.sig.txt breast 0 &> wgcna.bowtie.instance.BRCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.COADREAD.sig.txt colonrectum 0 &> wgcna.bowtie.instance.COADREAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.ESCASTAD.sig.txt esophagusstomach 0 &> wgcna.bowtie.instance.ESCASTAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.HNSC.sig.txt headneck 0 &> wgcna.bowtie.instance.HNSC.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.KICHKIRCKIRP.sig.txt kidney 0 &> wgcna.bowtie.instance.KICHKIRCKIRP.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.LIHC.sig.txt liver 0 &> wgcna.bowtie.instance.LICH.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.LUADLUSC.sig.txt lung 0 &> wgcna.bowtie.instance.LUAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.PRAD.sig.txt prostate 0 &> wgcna.bowtie.instance.PRAD.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.THCA.sig.txt thyroid 0 &> wgcna.bowtie.instance.THCA.log &
Rscript WGCNA.L1HS.instance.bytissue.r result.rsem.TET.bowtie.instance/ VSTcnt.UCEC.sig.txt endometrium 0 &> wgcna.bowtie.instance.UCEC.log &

##################### LAML ########################
Rscript L1HSstatistics.1tissue.r data/rsem data/TET.LAML &> TET.LAML.log 
Rscript gene2L1HS.r result.rsem.TET.LAML &> lmLAML.log 
