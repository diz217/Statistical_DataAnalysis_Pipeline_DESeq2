# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 11:45:54 2025

@author: dz33
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import subprocess
import sys
import glob

from Cleantitle_v3 import cleantitle
from DownloadCounts_v2 import Download_process

def norm_gse_id(raw):
    raw = raw.strip()
    digits = ''.join(c for c in raw if c.isdigit())
    if digits: return f'GSE{digits}'
    else: raise Exception(f'{raw} contains no digits') 

def run_r(script,*args):
    cmd = ["C:/Program Files/R/R-4.5.2/bin/Rscript.exe",script]+list(args)
    print('>>> run',' '.join(cmd))
    try:
        proc = subprocess.run(cmd,capture_output= True, text= True, check=True)
        if proc.stderr:
            print(proc.stderr)
    except subprocess.CalledProcessError as e:
        if e.stdout:
            print(e.stdout)
        if e.stderr:
            print(e.stderr)
            
def iter_run(gse_id,has_rep,multifile):
    if not multifile:
        print("Warning: Zero processed files.")
        sys.exit(0)
    for flpath in multifile:
        if has_rep:
            print(f"Analysis: DESeq2 on {flpath}")
            run_r("run_deseq2_v2.R",gse_id,flpath)
        elif '_val_' in flpath:
            print(f"Analysis: Limma on {flpath} (logFC and mean only)")
            run_r("run_limma_v2.R",gse_id,flpath)
        elif '_count_' in flpath:
            print(f"Analysis: Limma+Voom on {flpath} (logFC and mean only)")
            run_r("run_limmaVoom_v2.R",gse_id,flpath)
        else:
            print(f'Warning: {flpath} is unexpected. skip...')
            
def check_file(name,clauseA,clauseB):
    # name: str | list[str]
    if isinstance(name, list):
        files = []
        for pat in name:
            files.extend(glob.glob(pat))
    else:
        files = glob.glob(name)
    if not files:
        print(clauseA)
        return 1
    elif os.path.getsize(files[0])==0:
        print(clauseB)
        return 2
    else:
        return 0

if len(sys.argv) >=2:
    raw_id = sys.argv[1]
else:
    raw_id = input("请输入 GSE 编号：")

gse_id = norm_gse_id(raw_id)
print(f"GSE ID: {gse_id}")

run_r("save_gse_meta.R",gse_id)
code = check_file(f'{gse_id}_series*.txt','Series matrix download error!','0-Size series matrix !')
if code: sys.exit(code)
print(f'样本信息已下载: {gse_id}_series*.txt')

has_rep = cleantitle(gse_id)
code = check_file(f'{gse_id}_meta_gsm_condition.tsv','Meta Information Matrix not obtained !','Zero Experimental conditions !')
if code: sys.exit(code)
print(f'样本信息矩阵已处理，得到： {gse_id}_meta_gsm_condition.tsv')

dn_gse = Download_process(gse_id,has_rep,out_dir = './test/')
dn_gse.main_process()
check_list = [f'{dn_gse.out_dir}{gse_id}_count_gsm_matrix*.tsv',f'{dn_gse.out_dir}{gse_id}_val_gsm_matrix*.tsv']
code = check_file(check_list,'Count Matrix failed to process !','Zero counts !')
if code: sys.exit(code)
print(f'计数矩阵已下载，得到： {gse_id}_(count|val)_gsm_matrix*.tsv')

iter_run(gse_id,has_rep,dn_gse.filepaths)
print('统计分析已完成，位于： ./DESeq2/ or ./Limma/')
