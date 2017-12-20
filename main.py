import shrunken_centriod as sc
import matplotlib.pyplot as plt

parent = '/Users/zhongningchen/LungCancer/mi_RNA_Isoform_Data_nxp/miRNa'
file_miRNA_LUAD = '/LUAD_miRNA_nxp.csv'
file_miRNA_LUSC = '/LUSC_miRNA_nxp.csv'
file_list = [file_miRNA_LUAD, file_miRNA_LUSC]

miRNA = sc.shrunkenCentroid(parent, file_list)



def plot_meanDF(dmag):
    
    
    
    meanDF = miRNA.xikPrime(dmag*miRNA.d0)
    
    
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


def plot_meanDF_series(n, inc):
    
    for i in range(n):
        plot_meanDF((i + 1)*inc)
    
    return None



plot_meanDF_series(5, 25)
