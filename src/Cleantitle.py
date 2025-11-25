# -*- coding: utf-8 -*-
"""
Created on Sun Nov 23 21:45:10 2025

@author: Ding Zhang
"""

import pandas as pd
import os
import re
import sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def Cleantitle(gse_id):
    df = pd.read_csv(f'{gse_id}_series.txt',sep=r'\t',header=0)
    df_old = df.copy(deep=True)
    lis = []
    for item in df['title'].values:
        numstr = re.search('(\d+)',item).group(1)
        lis.append(int(numstr))
    unique = len(set(lis))==len(lis)
    if unique:
        cleanlis =[]
        for item in df['title'].values:
            t2 = re.sub('\d+','',item,count=1)
            cleanlis.append(t2)
        df['title'] = cleanlis
    
    df2 = pd.DataFrame()
    df2 = df.groupby('title').size().reset_index(name='count')
    if (df2['count']==1).any():
        sep = r"[,\s_\-./():;`]+"
        lis = []
        for item in df['title'].values:
            token0 = re.split(sep,item)
            tokens = [t for t in token0 if t!='']
            lis.append('.'.join(t for t in tokens[:-1]))  
        df['title'] = lis
        df2 = df.groupby('title').size().reset_index(name='count')
        if (df2['count']==1).any():
            print('no duplicates, check title')
            sys.exit(1)
    if 'description' in df.columns:
        df3 = df[["geo_accession","title","description"]]
    else:
        df3 = df[["geo_accession","title"]]
    df3['title0'] = df_old['title']
    conditions = set(df['title'])
    
    df3.to_csv(f"./{gse_id}_meta_gsm_condition.tsv",sep='\t',index=False)
    with open(f'{gse_id}_conditions.txt','w') as f:
        for c in sorted(conditions):
            f.write(c+'\n')