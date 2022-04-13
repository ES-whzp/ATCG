import csv, re, os
import argparse
import subprocess
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math

def makeDir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def runCMD(cmd):
    # print(cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE, encoding="utf-8",
                   timeout=10000)

def extractPAM(inputPath, outputPath):
    fE = open(inputPath, 'r')
    pattern = re.compile('%s([A-Z]+)%s' % (fiveFlankMotif, threeFlankMotif))
    n = 0
    PAMCountDic = {}
    for line in fE:
        if n % 4 == 1:
            try:
                s = pattern.search(line)
                PAM = s.group(1)
                PAMCountDic[PAM] = PAMCountDic.get(PAM, 0) + 1
            except:
                pass
        # if n > 10:
        #     break
        n += 1
    fE.close()
    print('%s %d reads have been saved to %s\n'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), n/4, conExtractPAMPath))

    fo = open(outputPath, 'w', newline='')
    foC = csv.writer(fo)
    for k, v in PAMCountDic.items():
        # print(k, v)
        foC.writerow([k, v, str(len(k))])
    fo.close()

def filternN(extractPAMPath, n):
    PAMdic = dict()
    sum = 0
    with open(extractPAMPath) as fhand:
        for line in fhand:
            lineS = line.split(',')
            if 'N' not in lineS[0] and int(len(lineS[0])) == 6:
                PAM = lineS[0][-n:]
                PAMdic[PAM] = PAMdic.get(PAM, 0) + int(lineS[1])
                sum += int(lineS[1])
    return(PAMdic, sum)

def dic2DF(dic, n, log):
    NSet = set()
    if n == 6:
        for a in ['A', 'T', 'C', 'G']:
            for b in ['A', 'T', 'C', 'G']:
                for c in ['A', 'T', 'C', 'G']:
                    NSet.add(a+b+c)
    elif n == 4:
        for a in ['A', 'T', 'C', 'G']:
            for b in ['A', 'T', 'C', 'G']:
                NSet.add(a + b)
    dfDic = dict()
    if log == 0:
        for k1 in NSet:
            dfDic[k1] = {}
            for k2 in NSet:
                PAM = k1+k2
                reads = dic.get(PAM, -1)
                dfDic[k1][k2] = int(reads)
    if log == 2:
        for k1 in NSet:
            dfDic[k1] = {}
            for k2 in NSet:
                PAM = k1+k2
                reads = dic.get(PAM, 1e-10)
                dfDic[k1][k2] = math.log2(float(reads))
    if log == 10:
        for k1 in NSet:
            dfDic[k1] = {}
            for k2 in NSet:
                PAM = k1+k2
                reads = dic.get(PAM, 1e-10)
                dfDic[k1][k2] = math.log10(float(reads))
    return dfDic

def plot(conDic, expDic, conSum, expSum, weblogoReadsPath, weblogoSavePath, foldChangeCSVPath, foldThreshold, heatmapSavePath, n, log):
    ratio = conSum / expSum
    readsLst = []
    foldChangeLst = []
    foldDic = {}

    print('%s filtering the reads according to foldThreshold_%d_%dN...'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), foldThreshold, n))
    for k, v in conDic.items():
        if k in expDic:
            fold = int(v)/(int(expDic[k])*ratio)
            readsLst.append([k, v, expDic[k], int(expDic[k])*ratio, fold])
            foldDic[k] = fold
            if fold >= foldThreshold:
                foldChangeLst.append([k, round(fold)])
        else:
            foldDic[k] = int(v) // ratio
            foldChangeLst.append([k, int(v) // ratio])
            readsLst.append([k, v, 0, 0, 'Depletion'])

    if not os.path.exists(foldChangeCSVPath):
        with open(foldChangeCSVPath, 'w') as f:
            fc = csv.writer(f)
            fc.writerow(['PAM', 'reads_con', 'reads_exp', 'reads_exp_normlize', 'fold_change'])
            fc.writerow(['All PAM', conSum, expSum, ratio, 'None'])
            for ele in readsLst:
                fc.writerow(ele)

    print('%s plotting the Weblogo according to foldThreshold_%d_%dN...'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), foldThreshold, n))
    if not os.path.exists(weblogoReadsPath):
        # os.remove(weblogoReadsPath)
        with open(weblogoReadsPath, 'a') as fhand:
            a = 0
            for element in foldChangeLst:
                a = a + element[1]
                for i in range(int(element[1])):
                    fhand.write(element[0]+'\n')
    if n == 6:
        annotate = '-6,-5,-4,-3,-2,-1'
    elif n == 4:
        annotate = '-4,-3,-2,-1'
    xlable = "\"5' PAM position\""
    plotWeblogoCMD = 'weblogo -f %s -o %s -F jpeg -A dna --xlabel %s --annotate %s --color-scheme classic --resolution 600'%\
                     (weblogoReadsPath, weblogoSavePath, xlable, annotate)
    runCMD(plotWeblogoCMD)

    plotWeblogoCMD = 'weblogo -f %s -o %s -F eps -A dna --xlabel %s --annotate %s --color-scheme classic  --resolution 600'%\
                     (weblogoReadsPath, weblogoSavePath.replace('.jpg', '.eps'), xlable, annotate)
    runCMD(plotWeblogoCMD)

    if not os.path.exists(heatmapSavePath):
        print('%s plotting the heatmap according to foldThreshold_%d_%dN_log_%s...'%(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), foldThreshold, n, log))
        dfConDic = dic2DF(conDic, n, log)
        dfExpDic = dic2DF(expDic, n, log)
        dfFoldDic = dic2DF(foldDic, n, log)
        allDic = {'Control':dfConDic, 'experiment':dfExpDic, 'fold change':dfFoldDic}
        heatmap(allDic, heatmapSavePath, n)

def heatmap(allDic, savePath, n):
    f, axs = plt.subplots(1, 3, figsize=(100, 20))
    plt.subplots_adjust(left=0.15, bottom=0.1, top=0.9, right=0.95, hspace=0.2, wspace=0.25)
    nLst = [0, 1, 2]
    x = 0
    # cmap = sns.diverging_palette(180, 10, as_cmap=True)
    # colors = {'A': '#feb2b1', 'C': '#14c7fe', 'T': '#14f485', 'G': '#f8ffa3'}
    for title, dic in allDic.items():
        df = pd.DataFrame(dic).sort_index(axis=1, ascending=True).sort_index().T
        ax = axs[nLst[x]]
        # print(title, df)
        if n == 6:
            ax.set_title(title, fontsize=24)
            ax.tick_params(labelsize=16, colors='black', labeltop=True, labelbottom=False)
        if n == 4:
            ax.set_title(title, fontsize=48)
            ax.tick_params(labelsize=32, colors='black', labeltop=True, labelbottom=False)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=-90)
        # xtable = ax.table(cellColours= '#f8ffa3',
        #                    cellLoc='center', loc='top')

        sns.heatmap(data=df, annot=False, fmt = "d", linewidths=0.25, linecolor="black", cmap='Blues',
                annot_kws={'size': 20, 'weight': 'bold', 'color': 'black'}, ax = ax)

        x += 1
    plt.savefig(savePath, dpi=300, bbox_inches='tight')
    # plt.show()

parser = argparse.ArgumentParser(description='PAM Depletion Analysis Pipeline. Including 3 main function: \n\t1. Extract the PAM sequence \n\t2. Filter the PAM sequence \n\t3. plot the heatmap and weblogo', formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-c', '--control', nargs=1, metavar='', help='fastq path', type=str)
parser.add_argument('-e', '--experiment', nargs=1, metavar='', help='fastq path', type=str)
parser.add_argument('-f', '--fFlank', nargs=1, metavar='', help="5'-flanking sequence of PAM", type=str)
parser.add_argument('-t', '--tFlank', nargs=1, metavar='', help="3'-flanking sequence of PAM", type=str)
parser.add_argument('-ft', '--foldTH', nargs=1, metavar='', help="Fold threshold for plotting weblogo, default is 10", default=10, type=int)
# parser.add_argument('-type', '--casType', help="casType determine the xlable of weblogo, e.g. 'cas9', 'cas12'", default=0 ,type=int)
# parser.add_argument('-n', '--nt', nargs=1, metavar='', help="nucleotide number in heatmap and weblogo, default is 6", default=6, type=int)
parser.add_argument("-o", '--output', nargs=1, metavar='', help="output directory path", type=str)

args = parser.parse_args()
conInputPath = args.control[0]
expInputPath = args.experiment[0]
# NN = args.n
# log = args.log
foldThreshold = args.foldTH[0]
saveDir = args.output[0]
fiveFlankMotif = args.fFlank[0]
threeFlankMotif = args.tFlank[0]
# casType = args.type

makeDir(saveDir)
conName = conInputPath.split('/')[-1]
expName = expInputPath.split('/')[-1]
conExtractPAMPath = os.path.join(saveDir, conName + '_ExtractPAM.csv')
expExtractPAMPath = os.path.join(saveDir, expName + '_ExtractPAM.csv')
if not os.path.exists(conExtractPAMPath):
    print('%s Extracting the control PAMs...'%time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    extractPAM(conInputPath, conExtractPAMPath)
if not os.path.exists(expExtractPAMPath):
    print('%s Extracting the experiment PAMs...'%time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    extractPAM(expInputPath, expExtractPAMPath)

for n in [4, 6]:
    conDic, conSum = filternN(conExtractPAMPath, n)
    expDic, expSum = filternN(expExtractPAMPath, n)
    foldChangeCSVPath = os.path.join(saveDir, 'con_%s_exp_%s_foldChange_%dN.csv' % (conName, expName, n))
    weblogoReadsPath = os.path.join(saveDir, 'con_%s_exp_%s_foldThreShold_%d_%dN_WebLogo.txt' % (conName, expName, foldThreshold, n))
    weblogoSavePath = os.path.join(saveDir, 'con_%s_exp_%s_foldThreShold_%d_%dN_WebLogo.jpg' % (conName, expName, foldThreshold, n))
    for log in [0, 10]:
        heatmapSavePath = os.path.join(saveDir, 'con_%s_exp_%s_%dN_log%d_heatmap.png' % (conName, expName, n, log))
        if not os.path.exists(heatmapSavePath):
            plot(conDic, expDic, conSum, expSum, weblogoReadsPath, weblogoSavePath, foldChangeCSVPath, foldThreshold, heatmapSavePath, n, log)
