#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd
import numpy as np
import os

pd.set_option('display.max_columns', 15)
pd.set_option('display.width', 200)
sns.set(context='paper')
sns.set_palette("colorblind")

#Load the original orthlog table
CommonDf = pd.read_csv('../CommonDf.tsv', sep='\t', index_col=0) #12621
#drop gene families with no overlap
CommonDf.dropna(inplace=True, subset=['overlap']) #12618
AaLs = [x for x in list(CommonDf.columns) if x.endswith('_aa')]


CommonDf['Min'] = CommonDf[AaLs].apply(lambda row: min(row), axis=1)
CommonDf['Max'] = CommonDf[AaLs].apply(lambda row: max(row), axis=1)

################# figure 2
from scipy.stats import linregress
linregress(CommonDf.Min, CommonDf.Max)

CommonDf['MaxDif'] = CommonDf[AaLs].apply(lambda row: (max(row) - min(row)), axis=1)
CommonDf['RelMaxDif'] = CommonDf['MaxDif']*100/CommonDf['Min']

CommonDf['RelOverlap'] = CommonDf['overlap']*100/CommonDf['Min']
CommonDf['%MinId'] = CommonDf['AbsId']*100/CommonDf['Min']

with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'axes.labelsize':14, 
                             'axes.titlesize': 18,
                             'xtick.labelsize': 12,
                             'ytick.labelsize': 12,}))):
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig, ((ax1, ax2, ax3, ax4)) = plt.subplots(nrows=1, ncols=4,figsize=(10,4))
    #figure 2a scatterplot for min and max distribution of each cluster
    sns.regplot(x="Min", y="Max", data=CommonDf, y_jitter=0.05, line_kws={'color':'red'}, ax=ax1)
    ax1.set(xticks=range(0,8001,2000), xlim=(0,8000), ylim=(0,8000), xlabel="Shortest Ortholog\n(aa residues)", ylabel="Longest Ortholog\n(aa residues)")
    plt.setp(ax1.get_xmajorticklabels(), rotation=30)
    ax1.text(-0.3, 1.2, 'a', transform=ax1.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right')
    #figure 2b %maxdiff histogram
    sns.histplot(CommonDf[CommonDf['RelMaxDif'] > 0]['RelMaxDif'], bins=range(0,15,1), kde=False, ax=ax2)
    #counts, bins, patches = ax3.hist(CommonDf[CommonDf['%MaxDif'] > 0]['%MaxDif'], bins=10)
    ax2.set(xticks=range(1,15,1), xlim=(0,9), xlabel="Length difference (%)", ylabel="Ortholog Families")
    plt.setp(ax2.get_xmajorticklabels(), rotation=30)
    ax2.text(-0.3, 1.2, 'b', transform=ax2.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    # figure 2c Rel-overlap
    sns.histplot(CommonDf['RelOverlap'], bins=range(90,101,1), kde=False, ax=ax3)
    ax3.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Aligment Saturation (%)", ylabel="Ortholog Families")
    plt.setp(ax3.get_xmajorticklabels(), rotation=30)
    ax3.text(-0.3, 1.2, 'c', transform=ax3.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    # figure 2d %Abs identiy
    sns.histplot(CommonDf['%AbsId'], bins=range(90,101,1), kde=False, ax=ax4)
    ax4.set(xticks=range(90,101,2), xlim=(90,100), xlabel="Idenity (%)\n(within aligned region)", ylabel="Ortholog Families")
    plt.setp(ax4.get_xmajorticklabels(), rotation=30)
    ax4.text(-0.3, 1.2, 'd', transform=ax4.transAxes,
      fontsize=30, fontweight='bold', va='top', ha='right') 
    plt.tight_layout(pad=1.08)
    plt.savefig('Fig2_OrthologFamilies.pdf')
    plt.savefig('Fig2_OrthologFamilies.png')

#Min in sp
for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Min']:
            CommonDf.at[Index, 'Min_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Min_sp'] = m
CommonDf.Min_sp.value_counts()

for Index, Row in CommonDf.iterrows():
    m = 0
    for Sp in AaLs:
        if Row[Sp] == Row['Max']:
            CommonDf.at[Index, 'Max_sp'] = Sp
            m += 1
    if m > 1: CommonDf.at[Index, 'Max_sp'] = m
CommonDf.Max_sp.value_counts()

################# figure 3 Multi sp substitution
CommonDf['MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))
CommonDf['%MultiSp_Subs'] = (CommonDf['Subs'] - CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs']].sum(axis=1))*100/CommonDf['overlap']

plt.close('all')
plt.clf()
plt.rcParams['pdf.fonttype'] = 'truetype'
fig, ax2 = plt.subplots(nrows=1, ncols=1,figsize=(3,5))
# figure 3g sum of all substitutions
CommonDf[['#1#_Subs', 'Convergent_Subs', 'OneInOutId', 'OnlyInGpId', 'OnlyOutGpId',
       'NoId',]].sum().plot.bar(ax=ax2)
ax2.set(ylabel="Total Number of sites", xticklabels=["#1#_Subs", 'Two identities - inconsistent', 'One identity - inconsistent', 'Only Ingroup identical', 'Only Outgroup identical',
       'No identity',])
plt.setp(ax2.get_xmajorticklabels(), rotation=30, ha='right')

plt.tight_layout(pad=1.08)
plt.savefig('Fig3g_MultiSpSubs.pdf')

################# figure 4 Departure from the average tree
with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'font.size': 14.0,
                             'axes.labelsize':14, 
                             'axes.titlesize': 40,
                             'xtick.labelsize': 14,
                             'ytick.labelsize': 14,
                             'legend.fontsize': 12,
                             }))):    
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig2 = plt.figure(figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
    #ax1 = fig2.add_subplot(spec2[0, 0])
    ax2 = fig2.add_subplot(spec2[0, 0])
    #1 RF metric
    PcDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf.tsv', sep='\t', index_col=0)
    BraScStatSigDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf/FdrBraScStatSigDf.tsv', sep='\t', index_col=0)
    BraScDf = pd.read_csv('../Subs_gblocks.aln/WithNull/PoissonCorrectedDf/BranchScoreDf.tsv', sep='\t', index_col=0)
    Insig = BraScDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull >= 0.05].index, :].TreeFull
    Sig = BraScDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull < 0.05].index, :].TreeFull

    #2 Histogram Tree length
    iSign = PcDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull < 0.05].index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).sort_values(ascending=False)
    iInSign = PcDf.loc[BraScStatSigDf[BraScStatSigDf.TreeFull >= 0.05].index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).sort_values(ascending=False)
    ax2.hist([iInSign, iSign], bins=np.linspace(0,0.07,71), stacked=True, label=['Insignificant', 'Significant'], color=['grey', 'black'])
    MyMean = PcDf.loc[BraScStatSigDf.index.intersection(PcDf.index), ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']].sum(axis=1).mean()
    ax2.axvline(MyMean, color='k',linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = ax2.get_ylim()
    ax2.text(MyMean*1.05, max_ylim*0.9, 'Mean:\n{:.3f}'.format(MyMean))
    ax2.set(xlabel='PC Tree Length', ylabel = 'Ortholog Families')
    ax2.legend(title='Deviation from the\nmean tree', loc='upper right', title_fontsize=12, fancybox=True, framealpha=0.5)
    # ax2.set_title('a', loc='left', fontweight='bold') 
    # ax2.text(-0.15, 1.1, 'b', transform=ax2.transAxes,
    #   fontsize=40, fontweight='bold', va='top', ha='right')
    plt.tight_layout(pad=1.08)
    plt.savefig('Fig4_TreeDeparture.pdf')
    plt.savefig('Fig4_TreeDeparture.png')

for Cols in PcDf.columns:
    CommonDf[Cols.replace('Subs', 'Pc')] = PcDf[Cols]


################# Supplemetal figure 1
CommonDf[CommonDf.HS_Subs == 0].index.to_series().to_csv('HS.null.txt', index=False, header=False)
CommonDf[CommonDf.PT_Subs == 0].index.to_series().to_csv('PT.null.txt', index=False, header=False)
CommonDf[CommonDf.NL_Subs == 0].index.to_series().to_csv('NL.null.txt', index=False, header=False)
CommonDf[CommonDf.GG_Subs == 0].index.to_series().to_csv('GG.null.txt', index=False, header=False)
CommonDf[CommonDf['#1#_Subs'] == 0].index.to_series().to_csv('P1.null.txt', index=False, header=False)

################# Supplemetal figure 2 Subs Count
with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'axes.labelsize':14, 
                             'axes.titlesize': 18,
                             'xtick.labelsize': 12,
                             'ytick.labelsize': 12,}))):    
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig, ax2 = plt.subplots(nrows=1, ncols=1,figsize=(6,5))
    sns.boxplot(data=CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', 'NL_Subs', '#1#_Subs',]], showfliers = False, ax=ax2)
    ax2.set(xlabel='', ylim=(0,11), ylabel = 'Substiutions', xticklabels=['Human', 'Chimpanzee', 'Gorilla','Gibbon', '#1#',])
    plt.tight_layout(pad=1.08)
    plt.savefig('Supp2_Subs.pdf')
    plt.savefig('Supp2_Subs.png')

    
    

################# figure 5 Branch-specific PC distance and Norm Distance
PcDf = pd.read_csv('../Subs_gblocks.aln/Recalculated/PoissonCorrectedDf.tsv', sep='\t', index_col=0)
#PC distance
with plt.rc_context(dict(sns.axes_style("white", rc = {
    'axes.spines.right': False,
    'axes.spines.top': False,}),
                         **sns.plotting_context("paper", rc={
                             'axes.labelsize':14, 
                             'axes.titlesize': 18,
                             'xtick.labelsize': 12,
                             'ytick.labelsize': 12,}))):    
    plt.close('all')
    plt.clf()
    plt.rcParams['pdf.fonttype'] = 'truetype'
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,figsize=(12,5))
    sns.boxplot(data=PcDf[['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs']], showfliers = False, ax=ax1)
    sns.stripplot(data=PcDf[['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs']], dodge = True, palette='bright', jitter=0.03, size =4, alpha=0.3, ax=ax1)
    ax1.set(xlabel='', ylabel = 'Recalculated PC Length', ylim= (0, 0.15), xticklabels=['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon', ])
    ax1.text(-0.15, 1.1, 'a', transform=ax1.transAxes,
      fontsize=40, fontweight='bold', va='top', ha='right') 
    NormDf =  pd.read_csv('../Subs_gblocks.aln/Recalculated/PoissonCorrectedDf/NormPcDf.tsv', sep='\t', index_col=0)
    sns.boxplot(data=NormDf[['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs']], showfliers = False, ax=ax2)
    sns.stripplot(data=NormDf[['HS_Subs', 'PT_Subs', '#1#_Subs', 'GG_Subs', 'NL_Subs']], dodge = True, palette='bright', jitter=0.03, size =4, alpha=0.3, ax=ax2)
    ax2.set(xlabel='', ylabel = 'Normalised Length', xticklabels=['Human', 'Chimpanzee','#1#', 'Gorilla','Gibbon',])
    ax2.text(-0.15, 1.1, 'b', transform=ax2.transAxes,
      fontsize=40, fontweight='bold', va='top', ha='right') 
    plt.tight_layout(pad=1.08)
    plt.savefig('Fig5_PCdistance.pdf')
    plt.savefig('Fig5_PCdistance.png')


#### Candidate and lists
BraScStatSigDf = pd.read_csv('../Subs_gblocks.aln/Recalculated/PoissonCorrectedDf/BraScStatSigDf.tsv', sep='\t', index_col=0)
NormPcDf = pd.read_csv('../Subs_gblocks.aln/Recalculated/PoissonCorrectedDf/NormPcDf.tsv', sep='\t', index_col=0)
NormSigDf = pd.read_csv('../Subs_gblocks.aln/Recalculated/PoissonCorrectedDf/NormScStatSigDf.tsv', sep='\t', index_col=0)

for Cols in NormPcDf.columns:
    CommonDf[Cols.replace('Subs', 'Norm')] = NormPcDf[Cols]

MeanBr = PcDf.mean()
MeanNorm = NormPcDf.mean()
HighAvgDf = pd.DataFrame()
Hominid4orMore = CommonDf[(CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs']] >= 4).any(axis=1)]
n=0
Cols = PcDf.columns
for Indx, row in PcDf.iterrows():
    HighAvgDf.at[Indx, 'Branches'] = ''
    for Col in Cols:
        if CommonDf.loc[Indx, ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs']].max() < 4: continue
        if BraScStatSigDf.at[Indx, Col.replace('_Subs', '')] >= 0.05: continue
        if NormSigDf.at[Indx, Col] >= 0.05: continue
        if row[Col] <=  MeanBr[Col]: continue
        if NormPcDf.at[Indx, Col] <= MeanNorm[Col]: continue
        HighAvgDf.at[Indx, 'Branches'] = f"{Col},{HighAvgDf.at[Indx, 'Branches']}"
        n += 1
HighAvgDf = HighAvgDf[HighAvgDf.Branches != '']
HighAvgDf.Branches.value_counts()

#Lower than average 
LowAvgDf = pd.DataFrame()
n=0
Cols = PcDf.columns
for Indx, row in PcDf.iterrows():
    LowAvgDf.at[Indx, 'Branches'] = ''
    for Col in Cols:
        if Col == '#1#_Subs': continue
        if CommonDf.loc[Indx, ['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs']].max() < 4: continue
        if BraScStatSigDf.at[Indx, Col.replace('_Subs', '')] >= 0.05: continue
        if NormSigDf.at[Indx, Col] >= 0.05: continue
        if row[Col] >=  MeanBr[Col]: continue
        if NormPcDf.at[Indx, Col] >= MeanNorm[Col]: continue
        LowAvgDf.at[Indx, 'Branches'] = f"{Col},{LowAvgDf.at[Indx, 'Branches']}"
        n += 1
LowAvgDf = LowAvgDf[LowAvgDf.Branches != '']
LowAvgDf.Branches.value_counts()


#One sp overlap
LowInOneSp = LowAvgDf[[x in ['HS_Subs,', 'PT_Subs,', 'GG_Subs,', 'NL_Subs,'] for x in LowAvgDf.Branches]]
HighInOneSp = HighAvgDf[[x in ['HS_Subs,', 'PT_Subs,', 'GG_Subs,', 'NL_Subs,'] for x in HighAvgDf.Branches]]
HighInOneSp.index.union(LowInOneSp.index)
HighInOneSp.index.intersection(LowInOneSp.index)
LowHighDf =  LowInOneSp.loc[HighInOneSp.index.intersection(LowInOneSp.index),:].copy()
LowHighDf = LowHighDf + HighInOneSp.loc[HighInOneSp.index.intersection(LowInOneSp.index),:]
LowHighDf.Branches.value_counts()

CommonDf['lower'] = LowInOneSp.replace(['HS_Subs,', 'PT_Subs,', 'GG_Subs,', 'NL_Subs,'], ['Human', 'Chimpanzee', 'Gorilla','Gibbon'])
CommonDf['higher'] = HighInOneSp.replace(['HS_Subs,', 'PT_Subs,', 'GG_Subs,', 'NL_Subs,'], ['Human', 'Chimpanzee', 'Gorilla','Gibbon'])
CommonDf.loc[HighAvgDf[HighAvgDf.Branches == '#1#_Subs,'].index, 'higher'] = '#1#'

# Subs in outliers
LowOrHighInOneSp = CommonDf.loc[HighInOneSp.index.union(LowInOneSp.index),['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs',]].copy()
LowOrHighInOneSp.sum()



CommonDf.loc[HighAvgDf[HighAvgDf.Branches == 'HS_Subs,'].index, ['Gene', 'Description', '%HS_Subs', 'HS_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%HS_Subs',ascending=False).to_csv('HighAvg_Human.tsv', sep='\t', header=['Gene', 'Description', 'Branch-specific % subs per site', 'Norm branch length','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighAvgDf[HighAvgDf.Branches == 'PT_Subs,'].index, ['PT', 'Gene', 'Description', '%PT_Subs', 'PT_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%PT_Subs',ascending=False).to_csv('HighAvg_Chimpanzee.tsv', sep='\t', header=['Chimpanzee transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Norm branch length','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighAvgDf[HighAvgDf.Branches == '#1#_Subs,'].index, ['Gene', 'Description', '%#1#_Subs', '#1#_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%#1#_Subs',ascending=False).to_csv('HighAvg_P1.tsv', sep='\t', header=['Gene', 'Description', 'Branch-specific % subs per site', 'Norm branch length','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighAvgDf[HighAvgDf.Branches == 'GG_Subs,'].index, ['GG', 'Gene', 'Description', '%GG_Subs', 'GG_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%GG_Subs',ascending=False).to_csv('HighAvg_Gorilla.tsv', sep='\t', header=['Gorilla transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Norm branch length','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[HighAvgDf[HighAvgDf.Branches == 'NL_Subs,'].index, ['NL', 'Gene', 'Description', '%NL_Subs', 'NL_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'lower']].sort_values(by= '%NL_Subs',ascending=False).to_csv('HighAvg_Gibbon.tsv', sep='\t', header=['Gibbon transcript', 'Gene', 'Description', 'Branch-specific % subs per site', 'Norm branch length','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Lower in'], index=True, index_label='Ortholog family / Human transcript')

MixedHigh = HighAvgDf[[x not in ['HS_Subs,', 'PT_Subs,', 'GG_Subs,', '#1#_Subs,', 'NL_Subs,'] for x in HighAvgDf.Branches]]
CommonDf.loc[MixedHigh.index, ['Gene', 'Description', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap']].to_csv('MixedHighAvg.tsv', sep='\t', header=['Gene', 'Description','Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat'], index=True, index_label='Ortholog family / Human transcript')

CommonDf.loc[LowAvgDf[LowAvgDf.Branches == 'HS_Subs,'].index, ['Gene', 'Description', '%HS_Subs', 'HS_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= 'HS_Norm',ascending=True).to_csv('LowAvg_Human.tsv', sep='\t', header=['Gene', 'Description', 'Human % subs per site', 'Human norm branch length', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowAvgDf[LowAvgDf.Branches == 'PT_Subs,'].index, ['PT', 'Gene', 'Description', '%PT_Subs', 'PT_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= 'PT_Norm',ascending=True).to_csv('LowAvg_Chimpanzee.tsv', sep='\t', header=['Chimpanzee transcript', 'Gene', 'Description', 'Chimp % subs per site', 'Chimp norm branch length', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowAvgDf[LowAvgDf.Branches == 'GG_Subs,'].index, ['GG', 'Gene', 'Description', '%GG_Subs', 'GG_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= 'GG_Norm',ascending=True).to_csv('LowAvg_Gorilla.tsv', sep='\t', header=['Gorilla transcript', 'Gene', 'Description', 'Gorilla % subs per site', 'Gorilla norm branch length', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
CommonDf.loc[LowAvgDf[LowAvgDf.Branches == 'NL_Subs,'].index, ['NL', 'Gene', 'Description', '%NL_Subs', 'NL_Norm', 'HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs', 'overlap','RelOverlap', 'higher']].sort_values(by= 'NL_Norm',ascending=True).to_csv('LowAvg_Gibbon.tsv', sep='\t', header=['Gibbon transcript', 'Gene', 'Description', 'Gibbon % subs per site', 'Gibbon norm branch length', 'Human subs', 'Chimp subs', 'Gorilla subs', '#1# subs', 'Gibbon subs', 'Align overlap', 'Align Sat', 'Higher in'], index=True, index_label='Ortholog family / Human transcript')
#CommonDf[(CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']] == 0).all(axis=1)].loc[:,['HS_gene_id', 'Gene', 'Description', 'overlap', 'Subs', '%AbsId', 'RelOverlap']].to_csv('NoSubs.tsv', sep='\t', header=['Human Gene Id', 'Gene', 'Description', 'Alignment overlap', 'Substitutions','% Identity', 'Align Sat'], index=True, index_label='Ortholog family / Human transcript')

#Remove null substituion rows
CommonDf = CommonDf[(CommonDf[['HS_Subs', 'PT_Subs', 'GG_Subs', '#1#_Subs', 'NL_Subs']] > 0).any(axis=1)]


#Manual curation
CommonDf.loc['ENST00000450565','Comment'] = 'Most divergent gene on the human branch'
CommonDf.loc['ENST00000301633','Comment'] = 'Gibbon misssing starting methionine'
CommonDf.loc['ENST00000391916','Comment'] = 'All substitutions were within a short block'
CommonDf.loc['ENST00000637878', 'Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000008938', 'Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000259845','Comment'] = 'Among the top 5 divergent genes on the human branch'
CommonDf.loc['ENST00000454136','Comment'] = 'Among the top 5 divergent genes on the human branch'

CommonDf.loc['ENST00000640237','Comment'] = 'First half of the alighment is unrelaible even after filtering'
CommonDf.loc['ENST00000228468','Comment'] = 'First half of the alighment is unrelaible even after filtering'
CommonDf.loc['ENST00000390654','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000265729','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000614943','Comment'] = 'Frame-shift'
CommonDf.loc['ENST00000375464','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000446510','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000373383','Comment'] = 'Many substitutions are within a short block'
CommonDf.loc['ENST00000372398','Comment'] = 'Many substitutions are within a short block'


CommonDf.loc['ENST00000542996','Comment'] = 'Most divergent gene on the chimpanzee branch'
CommonDf.loc['ENST00000296280', 'Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000380041', 'Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000342995','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'
CommonDf.loc['ENST00000343470','Comment'] = 'Among the top 5 divergent genes on the chimpanzee branch'

CommonDf.loc['ENST00000359741','Comment'] = 'Most divergent gene on the #1# branch'
CommonDf.loc['ENST00000397893','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000274520','Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000299191', 'Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000345165','Comment'] = 'Many gaps in the alignment'
CommonDf.loc['ENST00000568377','Comment'] = 'Many gaps in the alignment'
CommonDf.loc['ENST00000396124', 'Comment'] = 'Among the top 5 divergent genes on the #1# branch'
CommonDf.loc['ENST00000433976','Comment'] = 'Among the top 5 divergent genes on the #1# branch'

CommonDf.loc['ENST00000370177','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000375098','Comment'] = 'Short Alignment'
CommonDf.loc['ENST00000359878','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000342790','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000651546','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000649185','Comment'] = 'Many subsituions are within a short block'
CommonDf.loc['ENST00000263317','Comment'] = 'Many subsituions are within a short block'

CommonDf.loc['ENST00000255499','Comment'] = 'Most divergent gene on the gorilla branch'
CommonDf.loc['ENST00000432056', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000254976', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000613760', 'Comment'] = 'Among the top 5 divergent genes on the gorilla branch'
CommonDf.loc['ENST00000255977','Comment'] = 'Among the top 5 divergent genes on the gorilla branch'

CommonDf.loc['ENST00000513010','Comment'] = 'Most divergent gene on the gibbon branch'
CommonDf.loc['ENST00000345088', 'Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000398462', 'Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000393330','Comment'] = 'Among the top 5 divergent genes on the gibbon branch'
CommonDf.loc['ENST00000523047','Comment'] = 'Among the top 5 divergent genes on the gibbon branch'


#save the norm data frame
CommonDf.loc[PcDf.index, :].to_csv('../NormDf.tsv', sep='\t', header=True, index=False)









