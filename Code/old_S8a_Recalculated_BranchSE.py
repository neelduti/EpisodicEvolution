#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os, sys, pandas as pd, numpy as np
from statsmodels.stats import weightstats


CommonDf = pd.read_csv('../../CommonDf.tsv', sep='\t', 
                       usecols=['HS', 'overlap', 'HS_Subs', 'PT_Subs', 
                                'GG_Subs', 'NL_Subs', '#1#_Subs',], 
                       index_col = 'HS')
SubCols = [x for x in CommonDf.columns if x.endswith('_Subs')]

#Remove zero alignment families
CommonDf = CommonDf[CommonDf.overlap > 0]
#Remove zero branch-specific substitutions families
CommonDf = CommonDf[(CommonDf[SubCols] > 0).any(axis=1)]

#Add 1 subs on each branch
CommonDf[SubCols] = CommonDf[SubCols] + 1
#Add 5 to overlap length
CommonDf['overlap'] = CommonDf.overlap + 5

#Poisson correction after artifically adding one substituion on each branch
PcDf = CommonDf[SubCols].copy(deep=True)
PcDf = PcDf.apply(lambda row: -np.log(1-row/(CommonDf.at[row.name, 'overlap'])), axis=1)
PcDf.to_csv('PoissonCorrectedDf.tsv', sep='\t', header=True, index=True, index_label='Id') 

#Substitution position
SubsPosDf = pd.read_csv('../SubPosDf.tsv', sep='\t', index_col=0, )
SubsPosDf = SubsPosDf.loc[SubsPosDf.index.intersection(PcDf.index),:]
#Fill no substituion branches
SubsPosDf.fillna('-', inplace=True)

def SubsSer(SubsString):
    if SubsString == '-': return([])#empty list
    else:
        SubSer = [int(x) for x in SubsString.split(',')]
        return(SubSer)

StatSignDf = pd.DataFrame(dtype=np.float64)
BsDf = pd.DataFrame(dtype=np.float64)
NormBsDf = pd.DataFrame(dtype=np.float64)
NormHominidBsDf = pd.DataFrame(dtype=np.float64)

m=0
for Indx, row in SubsPosDf.iterrows():
    n=0
    SampleDf = pd.DataFrame()
    #Increase alignment length by five
    Length = int(row['Length']) + 5
    PosSer =  pd.Series(range(1, Length+1, 1))
    #Add one extra substituted position on each branch
    HSsubSer = pd.Series(SubsSer(row[ 'HS']) + [PosSer.iloc[-5]])
    PTsubSer = pd.Series(SubsSer(row[ 'PT']) + [PosSer.iloc[-4]])
    P1subSer = pd.Series(SubsSer(row[ '#1#'])+ [PosSer.iloc[-3]])
    GGsubSer = pd.Series(SubsSer(row[ 'GG']) + [PosSer.iloc[-2]])
    NLsubSer = pd.Series(SubsSer(row[ 'NL']) + [PosSer.iloc[-1]])
    while n<1000:
        
        SamplePosSer = PosSer.sample(n=Length, replace=True)
        PosValCount = SamplePosSer.value_counts()
        SampleDf.at[n, 'HS'] = PosValCount.loc[PosValCount.index.intersection(HSsubSer.values)].sum()
        SampleDf.at[n, 'PT'] = PosValCount.loc[PosValCount.index.intersection(PTsubSer.values)].sum()
        SampleDf.at[n, '#1#'] = PosValCount.loc[PosValCount.index.intersection(P1subSer.values)].sum()
        SampleDf.at[n, 'GG'] = PosValCount.loc[PosValCount.index.intersection(GGsubSer.values)].sum()
        SampleDf.at[n, 'NL'] = PosValCount.loc[PosValCount.index.intersection(NLsubSer.values)].sum()
        
        n += 1
    Branches = SampleDf.columns
    #Poisson correction
    SampleDf = -np.log(1-SampleDf/Length)
    NormDf = SampleDf.copy(deep=True)
    TreeSum = NormDf.sum(axis=1)
    NormDf = NormDf.apply(lambda row: row/TreeSum[row.name], axis=1)
    NormDf.replace(np.nan, 0)
    
    #standard error for each branch
    for Col in SampleDf.columns:
        SampleDf[f"Var_{Col}"] = (SampleDf[f"{Col}"] - SampleDf[f"{Col}"].mean())**2
        NormDf[f"Var_{Col}"] = (NormDf[f"{Col}"] - NormDf[f"{Col}"].mean())**2
        BsDf.loc[Indx, f"SE_{Col}"] = ((SampleDf[f"Var_{Col}"].sum())/999)**0.5
        NormBsDf.loc[Indx, f"SE_{Col}"] = ((NormDf[f"Var_{Col}"].sum())/999)**0.5
        if Col == 'NL':continue
    SampleDf["TreeLenFull"] = TreeSum
    SampleDf["Var_TreeLenFull"] = (SampleDf["TreeLenFull"] - SampleDf["TreeLenFull"].mean())**2
    BsDf.loc[Indx, "SE_TreeFull"] = ((SampleDf["Var_TreeLenFull"].sum())/999)**0.5
    
    
    m += 1
    # print(m)
    # if m > 5: break

BsDf.to_csv('BootStrapSE.tsv', sep='\t', header=True, index=True, index_label='Id')
NormBsDf.to_csv('BootStrapSE.Norm.tsv', sep='\t', header=True, index=True, index_label='Id')  
    
