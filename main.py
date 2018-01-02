import shrunken_centriod as sc
import matplotlib.pyplot as plt
import pandas as pd
import math



from sklearn.model_selection import train_test_split

def loadData2(parent, file_list):
    #loads all data into dataframe dfList
    #for miRNA Lung cancer, it loads LUAD and LUSC as dataframes each
    #and stores to a list called dfList

    dfList = []
    for file in file_list:
        df_temp = pd.read_csv(parent + file)
        
        
        temp_dict = {}
        for i in range(len(df_temp['Unnamed: 0'])):
            temp_dict[i] = df_temp['Unnamed: 0'][i]
        df_temp = df_temp.rename(index = temp_dict)
        
        del df_temp['Unnamed: 0']
        
        df_temp['Label'] = pd.Series()
        df_temp['Label'].fillna(file, inplace = True)
        
        dfList += [df_temp]

    DF = pd.concat(dfList)
    DF = DF.fillna(0)
    

    return DF


def loadData(parent, file_list):
    #loads all data into dataframe dfList
    #for miRNA Lung cancer, it loads LUAD and LUSC as dataframes each
    #and stores to a list called dfList

    dfList = []
    for file in file_list:
        df_temp = pd.read_csv(parent + file)
        del df_temp['Unnamed: 0']
        dfList += [df_temp]

    
    return dfList 

def plot_meanDF(miRNA, dmag):
    
    meanDF = miRNA.xikPrime(dmag*miRNA.d0)
    i = 0
    v = [None]*len(file_list)
    for i in range(len(file_list)):
        fig = plt.figure(figsize = (12,6))
        v[i] = fig.add_subplot(121)
        t = range(meanDF.shape[0])
        v[i].vlines(t,[0], meanDF[file_list[i]])
        v[i].set_xlabel('genes expressed')
        v[i].set_title(file_list[i])
    
    plt.show()


def plot_meanDF_series(n, inc):
    for i in range(n):
        plot_meanDF((i + 1)*inc)
    
    return None






def discriminant(miRNA, meanDF, xi_star):
    l0 = miRNA.dfList[0].shape[0]
    l1 = miRNA.dfList[1].shape[1]
    n = l0 + l1
    prior = {file_list[0]: l0/n, file_list[1]: l1/n}
    
    ugn_star = list(set(miRNA.ugn[-1]).intersection(xi_star.index[:]))
    temp_disc_list = []
    for c in file_list:
        disc_sum = 0
        for gene in ugn_star:        
            disc = (xi_star[gene] - meanDF[c][gene])**2/(miRNA.sList[gene] + miRNA.s_knot)**2 - math.log(prior[c])
            disc_sum += disc 
        temp_disc_list += [disc_sum]
    return file_list[temp_disc_list.index(min(temp_disc_list))]
    #return temp_disc_list

def splitList(dfList):

    trainList = []
    testList = []
    for df in dfList:
        train, test = train_test_split(df, test_size=0.1)
        trainList += [train]
        testList += [test]
    
    return trainList, testList

def crossValidation(fold = 10, delta = 1):
    
    accuracy = 0
    for f in range(fold):
        trainList, testList = splitList(dfList)
        miRNA = sc.shrunkenCentroid(trainList, file_list)
        meanDF = miRNA.xikPrime(delta*miRNA.d0)
        temp = 0
        n = 0
        
        for classes in range(len(testList)):
            
            for i in range(testList[classes].shape[0]):
                if discriminant(miRNA, meanDF, dfList[classes].iloc[i,]) == file_list[classes]:
                    temp += 1
                
            n += testList[classes].shape[0]
            
        temp = temp/n
        accuracy += temp
        print (temp)
    accuracy = accuracy/fold
    return accuracy

def oneFoldValidation(delta):
    trainList, testList = splitList(dfList)
    miRNA = sc.shrunkenCentroid(trainList, file_list)
    meanDF = miRNA.xikPrime(delta*miRNA.d0)
    temp = 0
    n = 0
    for classes in range(len(testList)):
        for i in range(testList[classes].shape[0]):
            if discriminant(miRNA, meanDF, dfList[classes].iloc[i,]) == file_list[classes]:
                temp += 1
            #print (discriminant(miRNA, meanDF, dfList[classes].iloc[i,]) == file_list[classes])
        n += testList[classes].shape[0]    
        
    return temp/n

parent = '/Users/zhongningchen/LungCancer/mi_RNA_Isoform_Data_nxp/miRNa'
file_miRNA_LUAD = '/LUAD_miRNA_nxp.csv'
file_miRNA_LUSC = '/LUSC_miRNA_nxp.csv'
file_list = [file_miRNA_LUAD, file_miRNA_LUSC]
dfList = loadData(parent, file_list)

i = 0
result = []
while i < 7:
    result += [i*5, crossValidation(fold = 10, delta = i*5)]
    print ('Done with ' + str(i*5))
    i += 1
