# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 10:13:39 2025

@author: dtrzL
"""
import requests
import re
import sys
import os
from pathlib import Path
import pandas as pd 
import gzip,tarfile,rarfile,zipfile,py7zr
import shutil
from scipy.io import mmread

class Download_process():
    def __init__(self,gse_id,has_rep,out_dir = './'):
        self.gse_id = gse_id
        self.rep = has_rep
        self.out_dir = out_dir
        os.makedirs(out_dir,exist_ok=True)
        self.int_set = {'raw','count','matrix','.mtx'}
        self.frag_set = {'rpkm','fpkm','tpm','cpm'}
        self.bad_set = {'ercc', 'spike', 'phix', 'rrna', 'gfp', 'rfp', 'mcherry', 'mt-', 'summary',
                        'report', 'statistics', 'sample_sheet', 'bed', 'bedgraph','peak','hepa', 'cluster'}
        self.comment_set = tuple(['!','#','notes','summary','%','^'])
        self.bad_ext_set = {'.bed','.bedgraph','.wig','bigwig','.bam','.sam','.gtf','.vcf','.soft',
                            '.gff','.peak','.narrowpeak','.bw','.fastq','.fq','.family','.idat','.sra',
                            '.ann','.rds','.rdata','.rda','.h5','.h5ad','.loom','.cel','.doc','.docx'}
        self.colname_which = None 
    def main_process(self,):
        num = int(self.gse_id[3:])
        block = num//1000
        prefix = f"GSE{block}nnn"
    
        base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{self.gse_id}/suppl/"
        resp = requests.get(base_url)
        if resp.status_code !=200:
            print(f"stauts: {resp.status_code}")
        html = resp.text
        hrefs = re.findall(r'href="([^"]+)"',html)
        files = [h for h in hrefs if not h.startswith('?') and not h.startswith("/") and not h.startswith("http")]
        print('Files are being inspected...')
        paths = []
        for fn in files:
            low = fn.lower()
            if any(s in low for s in self.bad_set):
                continue
            if any(s in low for s in self.frag_set) or any(s in low for s in self.int_set):  
                url = base_url+fn
                local_path = self.out_dir+fn
                if not Path(local_path).is_file():
                    r = requests.get(url,stream=True)
                    r.raise_for_status()
                    with open(local_path,'wb') as f:
                        for chunk in r.iter_content(chunk_size=1024*1024):
                            if chunk: f.write(chunk)
                paths.append(local_path)
        if not paths:
            print('No legit count files inspected by filenames. Pipelines cannot be continued. Please verify manually.')
            sys.exit(0)
        df_cond = pd.read_csv(f'{self.gse_id}_meta_gsm_condition.tsv',sep='\t',header=0)
        gsms = set(df_cond['geo_accession'])
        filepath_list = []
        for ii,path1 in enumerate(paths): 
            print(f'No.{ii}: {path1} in processing...')
            mul,flname = self.decompress_file(path1)
            if mul:
                cmb_df = self.gsm_cmb(gsms,flname) # assume 1st col is gene
            else:
                table = self.gsm_reader(str(flname),self.gse_id) # assume 1st col is gene
                if table is None or table.shape[0]==0 or table.shape[1]==0:
                    print(f'{str(flname)} does not return content. skip...')
                    continue
                table.columns = table.columns.astype(str).str.strip().str.strip('"').str.strip("'")
                df_cols = df_cond.drop(columns=['condition'])
                keep_cols,keep_gsm = [],[]
                for colname in table.columns:
                    mask = (df_cols==colname).any(axis=1)
                    if mask.any():
                        idx = df_cols.index[mask].to_list()
                        if len(idx)>1:
                            print(f'Warning: In {str(flname)}, column {colname} has multiple matches in the series matrix. Please pay attention. 1st match is used by default.')
                        keep_cols.append(colname); keep_gsm.append(df_cols.loc[idx[0],'geo_accession'])
                cmb_df = table[keep_cols]
                cmb_df.columns = keep_gsm
            if cmb_df is not None and not cmb_df.empty:
                # cmb_df gsm order == meta gsm order
                cmb_df = self.reorder(cmb_df, df_cond)
                cmb_df.index = cmb_df.index.str.strip().str.strip('"').str.strip("'")
                print('Count matrix successfully processed.')
                cmb_df = cmb_df.apply(pd.to_numeric,errors='coerce')
                has_nonnumeric = cmb_df.isna().any(axis=1)
                if sum(has_nonnumeric)>20:
                    print(f'Warning: {sum(has_nonnumeric)} gene rows have non-numeric values. This is unexpected and should be checked manually.')
                cmb_df_clean = cmb_df.loc[~has_nonnumeric,:]
                has_float = (cmb_df_clean.stack() %1 !=0).any()
                filename = f'{self.out_dir}{self.gse_id}_val_gsm_matrix_{ii}.tsv' if has_float else f'{self.out_dir}{self.gse_id}_count_gsm_matrix_{ii}.tsv'
                print(f'file: {filename}')
                cmb_df_clean = cmb_df_clean[~cmb_df_clean.index.duplicated(keep='first')]
                cmb_df_clean.to_csv(filename,sep='\t',index=True)
                filepath_list.append(filename)
        self.filepaths =  filepath_list if filepath_list else None 
    def decompress_file(self,path1):
        flname = Path(path1).with_suffix('') # rm the last suffix
        mul = 1
        if open(path1,'rb').read(2)==b'\x1f\x8b': #.gz
            mul = 0
            with gzip.open(path1,'rb') as f:
                with open(flname,'wb') as g:
                    shutil.copyfileobj(f,g)
        elif tarfile.is_tarfile(Path(path1)):
            with tarfile.open(path1, "r:*") as tar:
                tar.extractall(flname)
        elif rarfile.is_rarfile(Path(path1)):
            with rarfile.RarFile(path1) as rar:
                rar.extractall(flname)
        elif zipfile.is_zipfile(Path(path1)):
            with zipfile.ZipFile(path1) as z:
                z.extractall(flname)
        elif open(path1,'rb').read(6)==b'\x37\x7A\xBC\xAF\x27\x1C': #.7z
            with py7zr.SevenZipFile(path1,mode='r') as z:
                z.extractall(flname)
        else: # Not compressed
            mul=0
            flname = Path(path1)
        if mul:
            count = 0
            with os.scandir(flname) as cur: 
                for entry in cur: 
                    if entry.is_file() and not entry.name.startswith('.'):
                        count +=1
                    if count>1:
                        break
            mul = 0 if count==1 else 1
        return mul,flname
    
    def gsm_cmb(self,gsms,flname):
        gsm_pat = re.compile(r'(GSM\d+)',re.I)
        cmb_list = []
        with os.scandir(flname) as cur:
            for entry in cur: 
                if entry.is_file() and gsm_pat.search(entry.name):
                    gsm = gsm_pat.search(entry.name).group(1)
                    if gsm in gsms:
                        if any(entry.name.lower().endswith(c) for c in self.bad_ext_set):
                            continue
                        if any(c in entry.name.lower() for c in self.int_set) or any(c in entry.name.lower() for c in self.frag_set):
                            table = self.gsm_reader(entry.path,gsm) # dataframe
                            gsm_count_colInd = 0
                            if table is not None and table.shape[1]>1:
                                gsm_count_colInd = self.guess_count_col(table)
                            table = table.iloc[:,[gsm_count_colInd]]
                            table.columns = [gsm]
                            cmb_list.append(table)
        if cmb_list:
            cmb_df = pd.concat(cmb_list,axis=1)
            cmb_df.fillna(0,inplace=True)
            return cmb_df
        
    def gsm_reader(self,path1,gsm):
        with open(path1,'rb') as f:
            head1024 = f.read(1024)
        if b'\x00' in head1024: # is_binary   
            if head1024[:2]==b'\x1f\x8b': #.gz
                with gzip.open(path1,'rb') as f:
                    next_path = Path(path1).with_suffix('')
                    with open(next_path,'wb') as g:
                        shutil.copyfileobj(f,g)
                return self.gsm_reader(str(next_path),gsm)
            elif head1024[:2]==b'PK' and path1.endswith('.xlsx'):
                cnt = pd.read_excel(path1,header=None)
            elif head1024[:8] == b'\xD0\xCF\x11\xE0\xA1\xB1\x1A\xE1' and path1.endswith('.xls'):
                cnt = pd.read_excel(path1,header=None)
            else: 
                print(f'Warning: {path1} is an unexpected binary file, skip...')
                return None 
        elif head1024.startswith(b'%%MatrixMarket'): # matrix.mtx
            genes_cnt = pd.DataFrame()    
            with os.scandir(Path(path1).parent) as cur:
                for entry in cur:
                    if entry.is_file() and gsm.lower() in entry.name.lower() and 'feature' in entry.name.lower():
                        genes_cnt = self.gsm_reader(entry.path,gsm)
                        break
            if genes_cnt.shape==(0,0):
                print(f'Warning: genes not procured for {path1}')
                return None
            mtx = mmread(path1)
            is_integer = b'integer' in head1024
            sparse_mtx = mtx.tocsr()
            if sparse_mtx.shape[0]!=genes_cnt.shape[0]:
                print(f'Warning: features do not match with matrix rows for {path1}')
                return None 
            res = sparse_mtx.sum(axis=1) if is_integer else sparse_mtx.mean(axis=1) # np matrix
            cnt_cl = pd.DataFrame()
            cnt_cl[0],cnt_cl[1] = genes_cnt.index,res.A1
            return self.process_colname_index(cnt_cl)
        else: # utf-8
            if path1.endswith('.csv'):
                cnt = pd.read_csv(path1,header=None,sep=None,engine='python')
            elif path1.endswith('.tsv'):
                cnt = pd.read_csv(path1,header=None,sep=r'\t',engine='python')
                if not cnt.shape[1]>1:
                    cnt = pd.read_csv(path1,header=None,sep=None,engine='python')
            else: # txt, none-suffix, and user-defined suffix
                try:
                    cnt = pd.read_csv(path1,header=None,sep=None,quoting=3,engine='python')
                except Exception as e:
                    print(f'Warning: {path1} cannot be processed by pandas, exception: {e}')
                    return None
        nnz = cnt.notna().sum(axis=1)
        valid_nnz = nnz[nnz>0]
        if valid_nnz.empty:
            print(f'Warning: {path1} is either empty or full of nan values. skip...')
            return None
        len_col = valid_nnz.mode().iat[0]
        keep = (nnz==len_col)
        pos_1  = (keep==True).values.argmax()
        if pos_1>0:
            keep.iat[pos_1-1] = nnz.iat[pos_1-1]==len_col-1
        cnt_cl = cnt[keep].dropna(axis=1,how='all')
        cmt_mask = ~cnt_cl.iloc[:,0].astype(str).str.lower().str.startswith(self.comment_set)
        cnt_cl = cnt_cl[cmt_mask].reset_index(drop=True)
        return self.process_colname_index(cnt_cl)
    
    def process_colname_index(self,cnt_cl):
        na_count0 = pd.to_numeric(cnt_cl.iloc[0],errors='coerce').isna().sum()
        na_count1 = pd.to_numeric(cnt_cl.iloc[1],errors='coerce').isna().sum()
        if na_count0>na_count1:
            cnt_cl.columns = cnt_cl.iloc[0].values
            cnt_cl = cnt_cl.iloc[1:].reset_index(drop=True)
        cnt_cl.set_index(cnt_cl.columns[0],inplace=True)
        cnt_cl.index.name = 'Gene_ID'
        return cnt_cl
    
    def guess_count_col(self,table):
        self.count_col = {'count','matrix','rpkm','tpm','cpm','fpkm'}
        for i,colname in enumerate(table.columns):
            if any(c in colname.lower() for c in self.count_col):
                self.colname_which = i
                return i
        if not self.colname_which:
            input1 = input('Attention: please enter the column index of the counts in GSM files (0-indexed):')
            str_val = ''.join(s for s in input1 if s.isdigit())
            val = int(str_val) if str_val else None
            if not val or val>table.shape[1]:
                print(f'Warning: {val} exceeds the maximum column number. The last column will be used')
                val = table.shape[1]
            elif val==0:
                print('Unexpected values: 0. Please use the manual mode.')
                sys.exit(1)
            self.colname_which = val-1
            return val-1
        return self.colname_which
    def reorder(self,cmb_df,df_cond):
        if df_cond['geo_accession'].to_list()!=cmb_df.columns.to_list():
            len_cond,len_cmb = len(df_cond['geo_accession']),len(cmb_df.columns)
            if len_cmb<len_cond:
                df_cond = df_cond[df_cond['geo_accession'].isin(cmb_df.columns)]
            cmb_df = cmb_df.loc[:,df_cond['geo_accession'].tolist()]
            if len_cmb<len_cond:
                df_cond.to_csv(f"{self.gse_id}_meta_gsm_condition.tsv",sep='\t',index=False)
                df_ctrl = pd.read_csv(f'{self.gse_id}_control_experiments.tsv',sep='\t',engine='python',header=None)
                keep_cond = set(df_cond['condition'])
                df_ctrl = df_ctrl[df_ctrl[0].isin(keep_cond) & df_ctrl[1].isin(keep_cond)].reset_index(drop=True)
                df_ctrl.to_csv('{self.gse_id}_control_experiments.tsv',sep='\t',index=False)
        return cmb_df
        

