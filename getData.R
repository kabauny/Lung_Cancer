library(TCGAbiolinks)
library(DT)
library(DESeq)


query_isoform_LUSC <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Isoform Expression Quantification",
                  experimental.strategy = "miRNA-Seq"
                  )
GDCdownload(query_isoform_LUSC, method = "api", files.per.chunk = 10)
isoform_data_LUSC <- GDCprepare(query_isoform_LUSC)

query_isoform_LUAD <- GDCquery(project = "TCGA-LUAD",
                               data.category = "Transcriptome Profiling",
                               data.type = "Isoform Expression Quantification",
                               experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_isoform_LUAD, method = "api", files.per.chunk = 10)
isoform_data_LUAD <- GDCprepare(query_isoform_LUAD)



query_miRNA_LUSC <- GDCquery(project = "TCGA-LUSC",
                          data.category = "Transcriptome Profiling",
                          data.type = "Isoform Expression Quantification",
                          experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_miRNA_LUSC, method = "api", files.per.chunk = 10)
miRNA_data_LUSC <- GDCprepare(query_miRNA_LUSC)


query_miRNA_LUAD <- GDCquery(project = "TCGA-LUAD",
                             data.category = "Transcriptome Profiling",
                             data.type = "Isoform Expression Quantification",
                             experimental.strategy = "miRNA-Seq"
)
GDCdownload(query_miRNA_LUAD, method = "api", files.per.chunk = 10)
miRNA_data_LUAD <- GDCprepare(query_miRNA_LUAD)

write.csv(isoform_data_LUSC, 'isoform_data_LUSC.csv')
write.csv(isoform_data_LUAD, 'isoform_data_LUAD.csv')
write.csv(miRNA_data_LUSC, 'miRNA_data_LUSC.csv')
write.csv(miRNA_data_LUAD, 'miRNA_data_LUAD.csv')

