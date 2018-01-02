import pandas as pd
from sklearn import svm

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

parent = '/Users/zhongningchen/LungCancer/mi_RNA_Isoform_Data_nxp/miRNa'
file_miRNA_LUAD = '/LUAD_miRNA_nxp.csv'
file_miRNA_LUSC = '/LUSC_miRNA_nxp.csv'
file_list = [file_miRNA_LUAD, file_miRNA_LUSC]
combine = loadData2(parent, file_list)
Label = combine['Label']
Label.replace(file_list[0], 'LUAD', inplace = True)
Label.replace(file_list[1], 'LUSC', inplace = True)
X = combine.drop('Label', axis = 1)

