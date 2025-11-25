# -*- coding: utf-8 -*-outdir
"""
Created on Mon Nov 24 01:42:35 2025

@author: Ding Zhang
"""

import requests
import re
import sys
import os
from pathlib import Path
import pandas as pd 
import tarfile
import rarfile
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def download_process(gse_id,outdir='./'):
    num = int(gse_id[3:])
    block = num//1000
    prefix = f"GSE{block}nnn"

    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse_id}/suppl/"

    resp = requests.get(base_url)
    if resp.status_code !=200:
        print(f"stauts: {resp.status_code}")

    html = resp.text

    hrefs = re.findall(r'href="([^"]+)"',html)
    files = [h for h in hrefs if not h.startswith('?') and not h.startswith("/") and not h.startswith("http")]
    select = []
    for fn in files:
        low = fn.lower()
        if 'raw' in low or 'count' in low:
            select.append(fn)
    if not select:
        print('no raw or count files.')
        sys.exit(1)

    print('Files are being downloaded: ')
    for fn in select:
        print("",fn)

    paths = []
    for fn in select:
        url = base_url+fn
        local_path = outdir+fn
        r = requests.get(url,stream=True)
        r.raise_for_status()
        with open(local_path,'wb') as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                if chunk:
                    f.write(chunk)
        paths.append(local_path)
    for path1 in paths:
        df3_name = outdir+f"{gse_id}_meta_gsm_condition.tsv"
        df3 = pd.read_csv(df3_name,header=0,sep='\t',dtype=str)
        if is_singlefile(path1):
            counts = pd.read_csv(path1,sep=r'\s+',header=0,engine='python',compression='infer')
            counts.columns = counts.columns.str.strip().str.strip('"').str.strip('"')

            raw_sample_cols = list(counts.columns)

            
            if 'description' not in df3.columns:
                df3['description'] = df3['title']
            keep_cols = []
            keep2 = []
            for col in raw_sample_cols:
                if col in df3['geo_accession'].values:
                    gsm = col
                    keep_cols.append(col)
                    keep2.append(gsm)
                elif col in df3['description'].values:
                    gsm = df3.loc[df3['description']==col,'geo_accession'].iloc[0]
                    keep_cols.append(col)
                    keep2.append(gsm)
                elif col in df3['title0'].values:
                    gsm = df3.loc[df3['title0']==col,'geo_accession'].iloc[0]
                    keep_cols.append(col)
                    keep2.append(gsm)
            counts2 = counts[keep_cols]
            counts2.columns = keep2
            counts2.index = counts[counts.columns[0]]
            counts2.to_csv(f"./{gse_id}_count_gsm_matrix.tsv",sep='\t',index=True)
        else:
            Path(path1[:-3]).mkdir(parents=True,exist_ok=True)
            index = None
            if tarfile.is_tarfile(Path(path1)):
                with tarfile.open(path1, "r:*") as tar:
                    tar.extractall(path1[:-3])

            if rarfile.is_rarfile(Path(path1)):
                with rarfile.RarFile(path1) as rf:
                    rf.extractall(path1[:-3])

            counts2 = pd.DataFrame()
            for entry in os.scandir(path1[:-3]):
                if 'count' not in entry.name and 'raw' not in entry.name:
                    continue
                gsm = re.search(r'(gsm\d+)',entry.name,flags=re.IGNORECASE)
                if not gsm:
                    continue
                gsm1 = gsm.group(1)
                if gsm1 not in list(df3['geo_accession'].values):
                    continue
                small = pd.read_csv(f'{path1[:-3]}/{entry.name}',sep = r'[\s+,]',header=None, skiprows=[0],engine='python',dtype=str)
                counts2[gsm1] = small[small.columns[1]]
                index = small[small.columns[0]]
            counts2.index = index
            counts2.to_csv(f"./{gse_id}_count_gsm_matrix.tsv",sep='\t',index=True)
#        return counts2
            
        
def is_singlefile(path1):
    suf = Path(path1).suffixes
    if len(suf)>=1 and suf[-1]=='.gz' and '.tar' not in suf and '.rar' not in suf:
        return True
    else:
        return False 