#bash getrawdiscounts.sh L1HS
#bash getrawdiscounts.sh AluSx1
#bash getrawdiscounts.sh AluYa5
#bash getrawdiscounts.sh HERVK3-int
#bash getrawdiscounts.sh HERVK9-int
#bash getrawdiscounts.sh LTR5_Hs
# 
#Rscript plotdiscounts.r L1HS
#Rscript plotdiscounts.r AluSx1
#Rscript plotdiscounts.r AluYa5
#Rscript plotdiscounts.r HERVK3-int
#Rscript plotdiscounts.r HERVK9-int
#Rscript plotdiscounts.r LTR5_Hs
#
#find ../data/TET -name "*instance.cntTable" | grep -v v1  | grep -v discount > instancefiles
#parallel -j 20 -a instancefiles bash filterintronTET.sh
#rm instancefiles
#find ../data/TET -name "*cntTable2" | sed 's/\.\.\/data\/TET\///' | sort > cntTable2
#find ../data/TET -name "*discount.instance.cntTable" | grep -v v1 | sed 's/\.\.\/data\/TET\///' | sort > discountTable
#split -l 35 cntTable2  raw.
#split -l 35 discountTable  disc.
##rm cntTable2 discountTable
#
#Rscript discount.parallel.r raw.aa disc.aa &> aa.log &
#Rscript discount.parallel.r raw.ab disc.ab &> ab.log &
#Rscript discount.parallel.r raw.ac disc.ac &> ac.log &
#Rscript discount.parallel.r raw.ad disc.ad &> ad.log &
#Rscript discount.parallel.r raw.ae disc.ae &> ae.log &
#Rscript discount.parallel.r raw.af disc.af &> af.log &
#Rscript discount.parallel.r raw.ag disc.ag &> ag.log &
#Rscript discount.parallel.r raw.ah disc.ah &> ah.log &
#Rscript discount.parallel.r raw.ai disc.ai &> ai.log &
#Rscript discount.parallel.r raw.aj disc.aj &> aj.log &
#Rscript discount.parallel.r raw.ak disc.ak &> ak.log &
#Rscript discount.parallel.r raw.al disc.al &> al.log &
#Rscript discount.parallel.r raw.am disc.am &> am.log &
#Rscript discount.parallel.r raw.an disc.an &> an.log &
#Rscript discount.parallel.r raw.ao disc.ao &> ao.log &
#Rscript discount.parallel.r raw.ap disc.ap &> ap.log &
#Rscript discount.parallel.r raw.aq disc.aq &> aq.log &
#Rscript discount.parallel.r raw.ar disc.ar &> ar.log &
#Rscript discount.parallel.r raw.as disc.as &> as.log &
#Rscript discount.parallel.r raw.at disc.at &> at.log &
#
cat maxdiff.instance.raw.a*.txt > maxdiff.instance.txt &
cat L1HSdiff.instance.raw.a*.txt > L1HSdiff.instance.txt &
 

