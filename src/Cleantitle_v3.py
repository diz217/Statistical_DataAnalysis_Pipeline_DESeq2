# -*- coding: utf-8 -*-
"""
Created on Thu Dec 11 22:27:14 2025

@author: Ding Zhang
"""
import re
import pandas as pd
import os
import sys
from collections import defaultdict, Counter

def cleantitle(gse_id):
    series = check_series_name(gse_id)
    df = pd.read_csv(series,sep=r'\t',header=0,engine='python')
    df_cols = set(df.columns)
    missing = [c for c in ['title','geo_accession'] if c not in df_cols]
    if missing:
        raise Exception(f'Missing columns from the series matrix: {missing}. Abort...')
    titles0 = df['title'].to_list()
    unique = prelimiary_clean(df)
    titles = df['title'].to_list()
    sep = SepSymbol(titles)
    head = findExample(titles,sep) 
    print(f'Example title: {head}')
    if unique:
        deep_cleanI(head,titles,sep,1)
    if has_unique(titles):
        head = findExample(titles,sep) 
        deep_cleanII(head,titles,sep,1)
    if is_unique(titles):
        head = findExample(titles,sep) 
        deep_cleanI(head,titles,sep,0)
    if is_unique(titles):
        head = findExample(titles,sep) 
        deep_cleanII(head,titles,sep,0)
    cleanse_sep(titles,sep)
    df['condition'] = titles
    df['title'] = titles0
    if is_unique(titles):
        print('Warning: automatic cleansing of titles failed. Limma Pipeline will be attempted...')
        print('Generating control groups...')
        # conditions
        conditions = cond_gen(titles,1)
        with open(f'{gse_id}_conditions.tsv','w') as f:
            f.write('\n'.join(conditions))
        if len(conditions)>20:
            to_continue = input(f'Warning: {len(conditions)} unique conditions are detected, which might be abnormal. Do you wish to continue: ')
            if re.search(r'^(q(uit)?|n(o|ah)?|stop|end|cancel|0|false|exit|abort|off|)$',to_continue, re.I):
                print('Okay. System exit. Please use the prompt mode.')
                sys.exit(0)               
        # control experiments.     
        ctrls = smart_control(conditions,sep)
        if ctrls:
            print_df_cond(conditions,df,gse_id)
            #w = max(len(s) for pair in ctrls for s in pair)
            with open(f'{gse_id}_control_experiments.tsv','w') as f:
                for exp1,exp2 in ctrls:
                    f.write(f'{exp1}\t{exp2}\n')
            print('Control groups are auto-generated (unique keys).')
            return 0
        print('Only one condition is detected. Aborting..')
        sys.exit(1)
    else:
        print('Titles cleansed. Proceed for Deseq2 pipline...')
        # conditions
        conditions = cond_gen(titles,0)
        with open(f'{gse_id}_conditions.tsv','w') as f:
            f.write('\n'.join(conditions))
        # control experiments.     
        ctrls = smart_control(conditions,sep)
        if ctrls:
            print_df_cond(conditions,df,gse_id)
            #w = max(len(s) for pair in ctrls for s in pair)
            with open(f'{gse_id}_control_experiments.tsv','w') as f:
                for exp1,exp2 in ctrls:
                    f.write(f'{exp1}\t{exp2}\n')
            print('Control groups are auto-generated (with duplicated keys).')
            return 1
        else:
            conditions = cond_gen(titles,1)
            with open(f'{gse_id}_conditions.tsv','w') as f:
                f.write('\n'.join(conditions))
            ctrls = smart_control(conditions,sep)
            if ctrls:
                print_df_cond(conditions,df,gse_id)
                #w = max(len(s) for pair in ctrls for s in pair)
                with open(f'{gse_id}_control_experiments.tsv','w') as f:
                    for exp1,exp2 in ctrls:
                        f.write(f'{exp1}\t{exp2}\n')
                print('Control groups are auto-generated (unique+duplicated keys).')
                return 0
            print('Only one unique condition is detected. Aborting..')
            sys.exit(1)
def deep_cleanI(head,titles,sep,swi):
    sep_pat1 = re.compile(rf'(^|[{sep}#\s])(\d+)(?=([{sep}#\s]|$))')
    sep_list = sep_pat1.findall(head)
    if sep_list: print(f'Indices to be check: {[i for i in sep_list]}')
    sep_set = []
    for sel,_,ser in sep_list:
        if (sel,ser) not in sep_set:
            sep_set.append((sel,ser))
    sep_set.reverse()
    for sel,ser in sep_set:
        pat = re.compile(rf'{L(sel)}(\d+)(?={R(ser)})')
        dict1 = defaultdict(list)
        t_len_min = 99999
        for i,t in enumerate(titles):
            t_len = 0
            for m in pat.finditer(t):
                dict1[i].append([m.group(1),m.span(1)])
                t_len += 1
            t_len_min = min(t_len,t_len_min) if t_len!=0 else t_len_min
        CriteriaI(dict1,titles,t_len_min) if swi else CriteriaII(dict1,titles,t_len_min)
    return 

def deep_cleanII(head,titles,sep,swi):
    sep_pat1 = re.compile(rf'(^|[{sep}#\s])([a-zA-Z]+\d+)(?=([{sep}#\s]|$))')
    sep_list = sep_pat1.findall(head)
    if sep_list: print(f'Indices to be check: {[i for i in sep_list]}')
    sep_set = []
    for sel,hun,ser in sep_list:
        pfix = re.search(r'([a-zA-z]+)',hun).group(1)
        if (sel,pfix,ser) not in sep_set:
            sep_set.append((sel,pfix,ser))
    sep_set.reverse()
    for sel,pfix,ser in sep_set:
        pat = re.compile(rf'{L(sel)}{pfix}(\d+)(?={R(ser)})')
        dict1 = defaultdict(list)
        t_len_min = 99999
        for i,t in enumerate(titles):
            t_len = 0
            for m in pat.finditer(t):
                dict1[i].append([m.group(1),m.span(1)])
                t_len += 1
            t_len_min = min(t_len,t_len_min) if t_len!=0 else t_len_min 
        CriteriaI(dict1,titles,t_len_min) if swi else CriteriaII(dict1,titles,t_len_min)
    return 

def CriteriaI(dict1,titles,t_len_min):  # remove seq + rep 
    if len(dict1)==len(titles):
        count = t_len_min
        while t_len_min<99999 and count>0: # can only operate by t_len_min
            val_list = [int(val[count-1][0]) for _,val in dict1.items()]
            if len(set(val_list))==len(titles): # seq detect
                del_by_index(titles,dict1,count-1)
            elif sorted(list(set(val_list)))==list(range(min(val_list),min(val_list)+len(set(val_list)))): # rep detect
                del_by_index(titles,dict1,count-1)
            count -=1
            
def CriteriaII(dict1,titles,t_len_min): # remove index, aggressive last sort
    count = t_len_min
    while t_len_min<99999 and count>0:
        val_list = [int(val[count-1][0]) for _,val in dict1.items()]
        is_index = False
        for i,index in enumerate(val_list[:-3]): # 4 pit
            is_index = is_index or val_list[i:i+4]==list(range(index,index+4))
            if is_index:
                break
        if is_index:
            del_by_index(titles,dict1,count-1)
        count -=1
    return

def del_by_index(titles,dict1,index):
     for i,t in enumerate(titles):
         if i in dict1:
             spanL,spanR = dict1[i][index][1]
             titles[i] = t[:spanL]+t[spanR:]
     return
 
def L(sel):
    return '^' if sel== '' else re.escape(sel)

def R(ser):
    return '$' if ser== '' else re.escape(ser)       
            
def prelimiary_clean(df):
    #step 1: find rep
    pattern = re.compile(r'rep(licate)?[\s\t_,.:]*\d+',re.I)
    temp = []
    for val in df['title'].values:
        temp.append(pattern.sub('',val))
    if list(df['title']) != temp:
        print('Title: replicates detected in preliminary steps.')
        df['title'] = temp
    #step 2: rep or seq without prefix at val's end
    temp = []
    sep = r"[,\s_\-/();`]+#"
    for val in df['title'].values:
        list1 = re.split(sep,val.strip().strip('"').strip("'"))
        if list1 and list1[-1].isdigit():
            temp.append((list1[-1],'_'.join(t for t in list1[:-1] if t)))
    if temp and len(temp)==len(df):
        nums = sorted(list(set(int(x[0]) for x in temp)))
        if len(nums)==len(df): # seq
            print('Title: sequence number detected at the end of the titles.')
            df['title'] = [x[1] for x in temp]
        elif nums==list(range(nums[0],nums[0]+len(nums))): #rep
            print('Title: replicate number detected at the end of the titles.')
            df['title'] = [x[1] for x in temp]
    return has_unique(df['title'].to_list())

def SepSymbol(titles):
    sep_map = {}
    pattern = re.compile(r"[,\s_\-/();`#]" )
    for i,t in enumerate(titles):
        if i>10:
            break
        sep_list = pattern.findall(t)
        for sep in sep_list:
            sep_map[sep] = sep_map.get(sep,0)+1
    sep_tup = sorted(sep_map.items(),key=lambda x:x[1],reverse=True)
    return sep_tup[0][0]

def findExample(titles,sep): 
    res,max_len = None,0
    num_pat = re.compile(r'(\d+)')
    for t in titles:
        leng = len(num_pat.findall(t))
        if leng>max_len:
            max_len = leng
            res = t
    return res if res else titles[0]

def has_unique(titles):
    df2 = pd.DataFrame()
    df2[0] = titles
    df2_group = df2.groupby(0).size().reset_index(name='count') 
    return (df2_group['count']==1).any()

def is_unique(title):
    return len(set(title))==len(title)

def check_series_name(gse_id): # by default, the larger file is selected. 
    res = []
    file_pat = re.compile(rf'{gse_id}_series[_\d{1}]*.txt')
    for entry in os.scandir():
        if entry.is_file() and file_pat.search(entry.name):
            res.append(entry.name)
    if len(res)>1:
        print(f'{len(res)} series matrix detected: {res}')
        index = input('Which series matrix do you want to process (0-indexed): ')
        if re.search(r'^(q(uit)?|n(o)?|na|not applicable)$',index,re.I):
            sys.exit(0)
        else:
            index = ''.join(c for c in index if c.isdigit())
            if not index:
                raise Exception(f'"{index}" is not an integer')
            if int(index)>len(res)-1:
                raise Exception(f'"{index}" out of range: [0, {len(res)-1}]')
            return res[int(index)]
    elif not res:
        raise Exception('No series matrices detected! Make sure to download the series matrices first.')
    else:
        return res[0]
    
def cleanse_sep(titles,sep):
    for i,key in enumerate(titles):
        titles[i] = f'{sep}'.join(k.strip() for k in key.strip().split(sep) if k.strip())
    return

def cond_gen(titles,inclu):
   df2 = pd.DataFrame()
   df2[0] = titles
   if not inclu:
       df2_group = df2.groupby(0).size().reset_index(name='count').query('count>1').sort_values(by=0)
   else:
       df2_group = df2.groupby(0).size().reset_index(name='count').query('count>=1').sort_values(by=0)
   return set(df2_group[0])

def smart_control(conditions,sep):
    dict2 = {}
    ctrls = []
    for cond in conditions:
        conl = [k.strip() for k in cond.strip().split(sep) if k.strip()]
        dict2[tuple(conl)] = len(conl)
    sort_items = sorted(dict2.items(),key=lambda x:x[1])
    for i,(keyi,numi) in enumerate(sort_items):
        for j,(keyj,numj) in enumerate(sort_items[i+1:]):
            pair = sum((Counter(keyi) & Counter(keyj)).values())
            if pair>=numi-1:
                ctrls.append((sep.join(keyi),sep.join(keyj)))
    return ctrls

def print_df_cond(conditions,df,gse_id):
    df3 = df[df['condition'].isin(conditions)]
    cols = [c for c in ["geo_accession","title","condition","description"] if c in df.columns]
    df3 = df3[cols].copy()
    df3.to_csv(f"./{gse_id}_meta_gsm_condition.tsv",sep='\t',index=False)
