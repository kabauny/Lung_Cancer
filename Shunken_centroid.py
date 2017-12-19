#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 10:55:00 2017

@author: zhongningchen
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 11:12:24 2017

@author: kabauny
"""
import numpy as np
import pandas as pd

def centroid_class_list(dfList, gene):
    #input: dfList is a list of classes Ck in the form of lsit of dataframes
    #and gene of interest
    
    #output: 1. list of all centroid for gene i of each class k, and overall centroid
    #as the last element. 2. length of lists as above

    centroid_class = []
    centroid_class_length = []
    
    temp_sum = 0
    temp_length = 0
    
    for dl in dfList:
        s = sum(dl[gene])
        l = len(dl[gene])
        
        centroid_class += [s/l]
        centroid_class_length += [l]
        
        temp_sum += s
        temp_length += l

    centroid_class += [temp_sum/temp_length]
    centroid_class_length += [temp_length]
    
    return centroid_class, centroid_class_length

def within_class_STD(dfList, gene):
    
    xi_list, xi_length = centroid_class_list(dfList, gene)
    n = sum(xi_length)
    K_length = len(xi_list)
    
    s_prime = 0
    
    k = 0
    while k < K_length:
        temp_sum = 0
        for j in dfList[k][gene]:
            temp_sum += (j - xi_list[-1])**2
        s_prime += temp_sum
        k += 1

    s = 1/(n - K_length)*s_prime
    return np.sqrt(s)

def d_stat(dfList, gene, s0):
    dgc_List = []
    
    ccl, l = centroid_class_list(dfList, gene)
    s = within_class_STD(dfList, gene)
    i = 0
    while i < len(ccl) - 1:
        mk = np.sqrt(1/l[i] + 1/l[-1])
        dgc_List += [(ccl[i] - ccl[-1])/(mk*(s + s0))]
        i += 1
    return dgc_List


loc = '/Users/zhongningchen/LungCancer/mi_RNA_Isoform_Data_nxp/miRNa'
file_miRNA_LUAD = '/LUAD_miRNA_nxp.csv'
file_miRNA_LUSC = '/LUSC_miRNA_nxp.csv'

df_miRNA_LUAD = pd.read_csv(loc + file_miRNA_LUAD)
del df_miRNA_LUAD['Unnamed: 0']

df_miRNA_LUSC = pd.read_csv(loc + file_miRNA_LUSC)
del df_miRNA_LUSC['Unnamed: 0']

df_miRNA_List = [df_miRNA_LUAD, df_miRNA_LUSC]

centroid_class_list(df_miRNA_List, )
