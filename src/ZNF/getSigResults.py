'''
Get the LINE family/zinc finger gene combinations from each cancer type that are positively correlated with q-value < 0.000001
'''

import sys
import glob
import os.path

resultsfilepath = sys.argv[1]
outputfile = sys.argv[2]

sigresultsdict = {}

for file in glob.glob(resultsfilepath + '/*/*_positive_lm.txt'):
  canctype = os.path.dirname(file).split('/')[-1]
  linefam = os.path.basename(file).rstrip('_results_positive_lm.txt')
  if canctype not in sigresultsdict:
    sigresultsdict[canctype] = {}
  sigresultsdict[canctype][linefam] = {}
  with open(file, 'r') as results:
    header = results.readline()
    for line in results:
      linesp = line.rstrip('\n').split('\t')
      genename = linesp[0]
      genecoef = float(linesp[3])
      genepval = linesp[4]
      geneqval = float(linesp[5])
      if (geneqval < 1e-7 and abs(genecoef)>0.4):
        sigresultsdict[canctype][linefam][genename] = []
        sigresultsdict[canctype][linefam][genename].append(genecoef)
        sigresultsdict[canctype][linefam][genename].append(genepval)
        sigresultsdict[canctype][linefam][genename].append(geneqval)


with open(outputfile, 'w+') as output:
  output.write('LineFamily''\t''Gene''\t''CancType''\t''GeneCoef''\t''GenePval''\t''GeneQval')
  for canctype, linefams in sigresultsdict.iteritems():
    for linefam, genenames in sigresultsdict[canctype].iteritems():
      for genename, info in sigresultsdict[canctype][linefam].iteritems():
        if info:
          output.write('\n' + linefam + '\t' + genename + '\t' + canctype + '\t' + str(info[0]) + '\t' + info[1] + '\t' + str(info[2]))
