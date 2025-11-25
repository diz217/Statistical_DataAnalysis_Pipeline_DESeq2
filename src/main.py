# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 01:08:04 2025

@author: Ding Zhang
"""

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import subprocess
import sys

from Cleantitle import Cleantitle
from DownloadCount import download_process

def norm_gse_id(raw):
    raw = raw.strip().upper()
    if raw[:3] =='GSE':
        return raw
    digits = ''.join(c for c in raw if c.isdigit())
    return f'GSE{digits}'

def run_r(script,gse_id,outdir='./'):
    cmd = ["C:/Program Files/R/R-4.5.1/bin/Rscript.exe",script,gse_id,outdir]
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

def check_file(name,clauseA,clauseB):
    if not os.path.exists(name):
        print(clauseA)
        return 1
    elif os.path.getsize(name)==0:
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
code = check_file(f'{gse_id}_series.txt','Series matrix download error!','0-Size series matrix !')
if code: sys.exit(code)
print(f'样本信息已下载: {gse_id}_series.txt')

Cleantitle(gse_id)
code = check_file(f'{gse_id}_meta_gsm_condition.tsv','Meta Information Matrix not obtained !','Zero Experimental conditions !')
if code: sys.exit(code)
print(f'样本信息矩阵已处理，得到： {gse_id}_meta_gsm_condition.tsv')

download_process(gse_id)
code = check_file(f'{gse_id}_count_gsm_matrix.tsv','Count Matrix failed to download !','Zero counts !')
if code: sys.exit(code)
print(f'计数矩阵已下载，得到： {gse_id}_count_gsm_matrix.tsv')

run_r("run_deseq2.R",gse_id)
print('统计分析已完成，位于： ./DESeq2/')

