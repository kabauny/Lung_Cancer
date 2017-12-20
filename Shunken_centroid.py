#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 16:54:36 2017

@author: zhongningchen
"""

import numpy as np
import pandas as pd
import statistics as stat
import matplotlib.pyplot as plt

def reList(dfList, gene):
    short_List = []
    
    for dl in dfList:
        
        short_List += [dl[gene]]
    
    return short_List

def centroid_class(class_series):


    centroid_class = []

    temp_sum = 0
    temp_length = 0
    
    for cs in class_series:
        s = sum(cs)
        l = len(cs)
        
        centroid_class += [s/l]

        temp_sum += s
        temp_length += l

    centroid_class += [temp_sum/temp_length]

    
    return centroid_class

def overall_centroid_List(dfList, geneList):
    x_list = []
    for gene in geneList:
        newList = reList(dfList, gene)
        x_list += [centroid_class(newList)[-1]]
    
    return x_list

def within_class_STD(class_series):
    
    xi_list = centroid_class(class_series)
    
    s_prime = 0
    n = 0
    k = 0
    while k < len(xi_list) - 1:
        temp_sum = 0
        for xij in class_series[k]:
            temp_sum += (xij - xi_list[k])**2
        s_prime += temp_sum
        n += len(class_series[k])
        k += 1

    s = 1/(n - k)*s_prime
    return np.sqrt(s)

def s_List(dfList, geneList):
    sList =[]
    
    for g in geneList:
        class_series = reList(dfList, g)
        sList += [within_class_STD(class_series)]
    
    return sList

def mk_list(dfList, gene_default):
    mkList = []
    temp_list = []
    temp_sum = 0
    for dl in dfList:
        temp_list += [len(dl[gene_default])]
        temp_sum += len(dl[gene_default])
    
    for tl in temp_list:
        
        mkList += [np.sqrt(1/tl + 1/temp_sum)]
    
    
    return mkList 


def d_class(class_series, s0, mk_list):
    dClass = []
    
    xi_list = centroid_class(class_series)
    s = within_class_STD(class_series)
    i = 0
    while i < len(xi_list) - 1:

        dClass += [(xi_list[i] - xi_list[-1])/(mk_list[i]*(s + s0))]
        i += 1
    return dClass

def d_list(df, geneList, s0, mkList):
    
    dList = []
    
    for gene in geneList:
        class_series = reList(df, gene)
        dList += [d_class(class_series, s0, mkList)]
    dList = np.array(dList)
    dList = dList.transpose()
    dList = dList.tolist()
    return dList 

def d_prime(d, delta):

    d_temp = abs(d) - delta
    if d_temp <= 0:
        d_temp = 0

    return np.sign(d)*d_temp

#def xik_prime()

###############################################################################


    


def loadData(loc, file_List):
    #Set up my Data, List of miRNA is stored in df.miRNA_List    

    df_miRNA_LUAD = pd.read_csv(loc + file_miRNA_LUAD)
    del df_miRNA_LUAD['Unnamed: 0']
    
    df_miRNA_LUSC = pd.read_csv(loc + file_miRNA_LUSC)
    del df_miRNA_LUSC['Unnamed: 0']
    
    df_miRNA_List = [df_miRNA_LUAD, df_miRNA_LUSC]
    
    return df_miRNA_List
    

def uniqueGeneNames(df_list):
    UGN = []
    
    common = set(df_list[0]).intersection(df_list[1])
    UGN += [list(set(df_list[0])^common)]
    UGN += [list(set(df_list[1])^common)]
    UGN += [list(common)]
    
    return UGN


###############################################################################
parent = '/Users/zhongningchen/LungCancer/mi_RNA_Isoform_Data_nxp/miRNa'
file_miRNA_LUAD = '/LUAD_miRNA_nxp.csv'
file_miRNA_LUSC = '/LUSC_miRNA_nxp.csv'
file_list = [file_miRNA_LUAD, file_miRNA_LUSC]
df_miRNA_List = loadData(loc = parent, file_List = file_list)

def initialize():

    
    ugn = uniqueGeneNames(df_miRNA_List)
    
    xiList = overall_centroid_List(df_miRNA_List, ugn[-1])
    
    sList = s_List(df_miRNA_List, ugn[-1])
    s_knot = stat.median(sList)
    mkList = mk_list(df_miRNA_List, ugn[-1][42])
    
    dList = d_list(df_miRNA_List, ugn[-1], s_knot, mkList)
    d0 = stat.median(dList[0])
    d1 = stat.median(dList[1])

    return ugn, xiList, sList, s_knot, mkList, dList, d0, d1
    

def xikPrime(s_List, s0, mk_List, d_List, d):
    
    mean_df = pd.DataFrame()
    
    k = 0
    while k < len(mk_List):
        i = 0
        temp = []
        while i < len(s_List):
            dprime = d_prime(d_List[k][i], delta = d)
            temp += [mk_List[k]*(s_List[i] + s0)*dprime]
            #temp += [xi_List[i] + mk_List[k]*(s_List[i] + s0)*dprime]
            i += 1
        mean_df[file_list[k]] = pd.Series(temp)
        k += 1
    
    return mean_df





def plot_meanDF(dmag):
    
    
    
    meanDF = xikPrime(s_List = sList, s0 = s_knot, mk_List = mkList, d_List = dList, d = d0*dmag)
    
    
    fig = plt.figure(figsize = (12,6))
    v_LUAD = fig.add_subplot(121)
    t = range(meanDF.shape[0])
    v_LUAD.vlines(t,[0], meanDF[file_list[0]])
    v_LUAD.set_xlabel('genes expressed')
    v_LUAD.set_title('LUAD miRNA Expression')
    
    fig = plt.figure(figsize = (12,6))
    v_LUSC = fig.add_subplot(121)
    t = range(meanDF.shape[0])
    v_LUSC.vlines(t,[0], meanDF[file_list[1]])
    v_LUSC.set_xlabel('genes expressed')
    v_LUSC.set_title('LUSC miRNA Expression')
    
    plt.show()
    
    plt.savefig('LUAD_LUSC_miRNA_' + str(d0*dmag) + '.png')


def plot_meanDF_series(n, inc):
    
    for i in range(n):
        plot_meanDF((i + 1)*inc)
    
    return None

ugn, xiList, sList, s_knot, mkList, dList, d0, d1 = initialize()

plot_meanDF_series(5, 25)



