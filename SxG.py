#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 21:25:57 2017

@author: zhongningchen
"""
import pandas as pd 

file_loc = '/Users/zhongningchen/LungCancer/miRNA_Isoform_Data/'

###############################################################################

file_isoform_LUAD = 'isoform_data_LUAD.csv'
df_isoform_LUAD = pd.read_csv(file_loc + file_isoform_LUAD)
iso_LUAD = df_isoform_LUAD[['isoform_coords', 'read_count', 'barcode']]


###############################################################################

UGi = list(set(iso_LUAD['isoform_coords']))
USi = list(set(iso_LUAD['barcode']))

###############################################################################
def createEmptyDF(geneList, barcodeList):
    df = pd.DataFrame()
    
    i = 0
    while i < len(geneList):
        df[geneList[i]] = pd.Series([None]*len(barcodeList))
        i += 1
        if i %500 == 0:
            print (i/len(geneList))
    print ('done with creating the empty array with columns name')
    i = 0
    temp_dict = {}
    while i < len(barcodeList):
        temp_dict[i] = barcodeList[i]        
        i += 1
        if i %500 == 0:
            print (i/len(barcodeList))
    print ('done creating barcode dictionary')
    df.rename(index = temp_dict, inplace = True)
    return df


'''
df = createEmptyDF(geneList =  UGi, barcodeList = USi)
df.to_csv('empty_df_LUAD.csv')
'''

###############################################################################

df = pd.read_csv('empty_df_LUAD.csv')


temp_dict = {}
i = 0
while i < len(USi):
    temp_dict[i] = USi[i]        
    i += 1

df.rename(index = temp_dict, inplace = True)
del df['barcode'] #note, it will be Unnamed: 0 rather than barcode in other cases


for s in USi:
    temp = iso_LUAD.loc[iso_LUAD['barcode'] == s]
    temp_index = temp.index.values
    
    for ti in temp_index:
        
        df.loc[s, temp.loc[ti, 'isoform_coords']] = temp.loc[ti, 'read_count']
    
    print ('done with ' + s)

df.fillna(0, inplace = True)
df.to_csv('LUAD_isoform_nxp.csv')

###############################################################################

###############################################################################


###############################################################################

file_isoform_LUSC = 'isoform_data_LUSC.csv'
df_isoform_LUSC = pd.read_csv(file_loc + file_isoform_LUSC)
iso_LUSC = df_isoform_LUSC[['isoform_coords', 'read_count', 'barcode']]


###############################################################################

UGi = list(set(iso_LUSC['isoform_coords']))
USi = list(set(iso_LUSC['barcode']))

###############################################################################



df = createEmptyDF(geneList =  UGi, barcodeList = USi)
df.to_csv('empty_df_LUSC.csv')


###############################################################################

df = pd.read_csv('empty_df_LUSC.csv')


temp_dict = {}
i = 0
while i < len(USi):
    temp_dict[i] = USi[i]        
    i += 1

df.rename(index = temp_dict, inplace = True)
del df['Unnamed: 0'] #note, it will be Unnamed: 0 rather than barcode in other cases


for s in USi:
    temp = iso_LUAD.loc[iso_LUAD['barcode'] == s]
    temp_index = temp.index.values
    
    for ti in temp_index:
        
        df.loc[s, temp.loc[ti, 'isoform_coords']] = temp.loc[ti, 'read_count']
    
    print ('done with ' + s)

df.fillna(0, inplace = True)
df.to_csv('LUSC_isoform_nxp.csv')

###############################################################################
###############################################################################
###############################################################################


file_miRNA_LUAD = 'miRNA_data_LUAD.csv'
df_miRNA_LUAD = pd.read_csv(file_loc + file_miRNA_LUAD)
miRNA_LUAD = df_miRNA_LUAD[['miRNA_ID', 'read_count', 'barcode']]


###############################################################################

UGi = list(set(miRNA_LUAD['miRNA_ID']))
USi = list(set(miRNA_LUAD['barcode']))

###############################################################################



df = createEmptyDF(geneList =  UGi, barcodeList = USi)
df.to_csv('empty_df.csv')


###############################################################################

df = pd.read_csv('empty_df.csv')


temp_dict = {}
i = 0
while i < len(USi):
    temp_dict[i] = USi[i]        
    i += 1

df.rename(index = temp_dict, inplace = True)
del df['Unnamed: 0'] 


for s in USi:
    temp = miRNA_LUAD.loc[miRNA_LUAD['barcode'] == s]
    temp_index = temp.index.values
    
    for ti in temp_index:
        
        df.loc[s, temp.loc[ti, 'miRNA_ID']] = temp.loc[ti, 'read_count']
    
    print ('done with ' + s)

df.fillna(0, inplace = True)
df.to_csv('LUAD_miRNA_nxp.csv')

###############################################################################




###############################################################################
###############################################################################


file_miRNA_LUSC = 'miRNA_data_LUSC.csv'
df_miRNA_LUSC = pd.read_csv(file_loc + file_miRNA_LUSC)
miRNA_LUSC = df_miRNA_LUSC[['miRNA_ID', 'read_count', 'barcode']]


###############################################################################

UGi = list(set(miRNA_LUSC['miRNA_ID']))
USi = list(set(miRNA_LUSC['barcode']))

###############################################################################



df = createEmptyDF(geneList =  UGi, barcodeList = USi)
df.to_csv('empty_df.csv')


###############################################################################

df = pd.read_csv('empty_df.csv')


temp_dict = {}
i = 0
while i < len(USi):
    temp_dict[i] = USi[i]        
    i += 1

df.rename(index = temp_dict, inplace = True)
del df['Unnamed: 0'] #note, it will be Unnamed: 0 rather than barcode in other cases


for s in USi:
    temp = miRNA_LUSC.loc[miRNA_LUSC['barcode'] == s]
    temp_index = temp.index.values
    
    for ti in temp_index:
        
        df.loc[s, temp.loc[ti, 'miRNA_ID']] = temp.loc[ti, 'read_count']
    


df.fillna(0, inplace = True)
df.to_csv('LUSC_miRNA_nxp.csv')

###############################################################################







