#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 16:54:36 2017

@author: zhongningchen
"""

import numpy as np
import pandas as pd
import statistics as stat



class shrunkenCentroid:
    
    
    def __init__(self, parent, file_list):
        #input: parent directory and name of files for analysis
        #initializing all necessary variables, see each function for details
        
        self.parent = parent
        self.file_list = file_list
        
        self.loadData() 
        self.uniqueGeneNames() 
        self.overall_centroid_List()
        self.s_List()
        self.s_knot = stat.median(self.sList)
        self.mk_list()
        
        self.d_list()
        self.d0 = stat.median(self.dList[0])
        #self.d1 = stat.median(self.dList[1]) 
        

    def loadData(self):
        #loads all data into dataframe dfList
        #for miRNA Lung cancer, it loads LUAD and LUSC as dataframes each
        #and stores to a list called dfList

        self.dfList = []
        for file in self.file_list:
            df_temp = pd.read_csv(self.parent + file)
            del df_temp['Unnamed: 0']
            self.dfList += [df_temp]
    
    
    def uniqueGeneNames(self):
        self.ugn = []
        #finds unique gene names between two data set
        #for miRNA Lung cancer, LUAD and LUSC has different lists of miRNA
        # that are expressed. This function finds all the common miRNA's 
        # and the stores as the last element of ugn. 
        # first two elements are miRNA unique LUAD and LUSC respectively
        common = set(self.dfList[0]).intersection(self.dfList[1])
        self.ugn += [list(set(self.dfList[0])^common)]
        self.ugn += [list(set(self.dfList[1])^common)]
        self.ugn += [list(common)]
        
  
    
    def reList(self, gene):
        short_List = []
        #given gene of interest, it returns the Series only for that gene 
        # for each calss. minimizes memory usage for computations below
        
        for dl in self.dfList:
            
            short_List += [dl[gene]]
        
        return short_List
    
    def centroid_class(self, class_series):
        #in: given Series of genes of interest for each class, 
        #retuns: the centroid of that gene with given class (xik), as well as 
        #xi as the last element of the list. 
    
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
    
    
    
    def overall_centroid_List(self):
        #from the dataframe of size |ugn| x |file_list|
        # it finds xi for all genes. and stores
        # in global variable xiList size |ugn| 
        self.xiList = []
        for gene in self.ugn[-1]:
            class_series = self.reList(gene)
            self.xiList += [self.centroid_class(class_series)[-1]]
        

    
    def within_class_STD(self, class_series):
        #given the Series for a specific gene size |ugn| x |file_list|
        # returns the within-class s, float
        xi_list = self.centroid_class(class_series)
        
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
    
    def s_List(self):
        #generates the within-class s for all genes 
        #and stores as global variable sList, size |ugn| 
        self.sList =[]
        for g in self.ugn[-1]:
            class_series = self.reList(g)
            self.sList += [self.within_class_STD(class_series)]
        

    
    def mk_list(self):
        #computes the mk values for the two classes and stores as list in global
        #variable mkList size |file_list|
        self.mkList = []
        temp_list = []
        temp_sum = 0
        for dl in self.dfList:
            temp_list += [dl.shape[0]]
            temp_sum += dl.shape[0]
        
        for tl in temp_list:
            
            self.mkList += [np.sqrt(1/tl + 1/temp_sum)]
    


    def d_class(self, class_series):
        #input: series of specific gene. size |ugn| x |file_list|
        #return the standarization value dik for specific gene, and class 
        # size 1 x |file_list|
        dClass = []
        
        xi_list = self.centroid_class(class_series)
        s = self.within_class_STD(class_series)
        i = 0
        while i < len(xi_list) - 1:
            dClass += [(xi_list[i] - xi_list[-1])/(self.mkList[i]*(s + self.s_knot))]
            i += 1
        
        return dClass
    
    def d_list(self):
        #generates all dik for all genes as dList size = |ugn| x |file_list|
        self.dList = []
        
        for gene in self.ugn[-1]:
            class_series = self.reList(gene)
            self.dList += [self.d_class(class_series)]
     
        
        self.dList = np.array(self.dList)
        self.dList = self.dList.transpose()
        self.dList = self.dList.tolist()
        
    

    def d_prime(self, d, delta):
        #input: d, float,  standarization value 
        #       delta, float, change in d 
        #output: float of same sign as d decreased by magnitude of delta
        d_temp = abs(d) - delta
        if d_temp <= 0:
            d_temp = 0
    
        return np.sign(d)*d_temp

    
    def xikPrime(self, Delta):
        #input: Delta, float, delta as above
        #output: df of size |ugn| x |file_list|, computes the modified mean 
        mean_df = pd.DataFrame()
        
        k = 0
        while k < len(self.mkList):
            i = 0
            temp = []
            while i < len(self.sList):
                dprime = self.d_prime(self.dList[k][i], delta = Delta)
                temp += [self.mkList[k]*(self.sList[i] + self.s_knot)*dprime]
                #temp += [self.xiList[i] + self.mkList[k]*(self.sList[i] + self.s_knot)*dprime]
                i += 1
            mean_df[self.file_list[k]] = pd.Series(temp)
            k += 1
        
        return mean_df





