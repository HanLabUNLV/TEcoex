'''
Get the LINE instances/zinc finger gene combinations from each cancer type that are positively correlated with q-value < 1e-7
'''

import sys
import glob
import os.path

resultsfilepath = sys.argv[1]
outputfile = sys.argv[2]

sigresultsdict = {}

for file in glob.glob(resultsfilepath + '*/*_positive_lm.txt'):
  canctype = os.path.dirname(file).split('/')[-1]
  gene_fam = os.path.basename(file).rstrip('_results_positive_lm.txt')
  (zincgene, linefam) = gene_fam.split('_', 1)
  if canctype not in sigresultsdict:
    sigresultsdict[canctype] = {}
  if zincgene not in sigresultsdict[canctype]:
    sigresultsdict[canctype][zincgene] = {}
  with open(file, 'r') as results:
    header = results.readline()
    for line in results:
      linesp = line.rstrip('\n').split('\t')
      linelocus = linesp[0]+":"+linefam
      genecoef = linesp[3]
      genepval = linesp[4]
      geneqval = float(linesp[5])
      if (linesp[6] == "NA"):
        continue;
      geneeta2 = float(linesp[6])
      if (geneqval < 1e-7 and abs(geneeta2)>0.4):
        sigresultsdict[canctype][zincgene][linelocus] = []
        sigresultsdict[canctype][zincgene][linelocus].append(genecoef)
        sigresultsdict[canctype][zincgene][linelocus].append(genepval)
        sigresultsdict[canctype][zincgene][linelocus].append(geneqval)
        print (canctype, zincgene, linelocus, geneqval)
        

print('Writing results')
with open(outputfile, 'w+') as output:
  output.write('LineFamily''\t''Gene''\t''CancType''\t''GeneCoef''\t''GenePval''\t''GeneQval')
  for canctype, zincgenes in sigresultsdict.iteritems():
    for zincgene, linefams in sigresultsdict[canctype].iteritems():
      for linefam, info in sigresultsdict[canctype][zincgene].iteritems():
        if info:
          output.write('\n' + linefam + '\t' + zincgene + '\t' + canctype + '\t' + info[0] + '\t' + info[1] + '\t' + str(info[2]))
