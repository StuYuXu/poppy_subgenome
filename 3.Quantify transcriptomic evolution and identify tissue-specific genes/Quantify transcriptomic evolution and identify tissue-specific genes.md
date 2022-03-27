# 1. Quantify transcriptomic evolution and identify tissue-specific genes

## 1.1 Single-copy orthologous genes

```shell
orthofinder -f three_Papaver_species_protein -S diamond -M msa  -t 48
```



## 1.2 Duplicates

```
Note, the expression data of P. somniferum and P. rhoeas were obtained from:
Yang, X., Gao, S., Guo, L. et al. Three chromosome-scale Papaver genomes reveal punctuated patchwork evolution of the morphinan and noscapine biosynthesis pathway. Nat Commun 12, 6030 (2021). https://doi.org/10.1038/s41467-021-26330-8
```



### 1.2.1 Normalized expression data

```R
## HN1 means P. somniferum
## YMR means P. rhoeas

HN1 <- read.csv("/RNASEQ/HN1.matrix",sep="\t",
                stringsAsFactors = F,header = T)

HN1_1 <- as.matrix(HN1[,-1])

library(preprocessCore)
HN1_1<-normalize.quantiles(HN1_1)

HN1[,2:7] <- HN1_1[,1:6]

m.max = apply(HN1[,2:7],1,max)
quantile(m.var, probs = seq(0, 1, 0.01))
HN1 = HN1[which(m.max>=1),]


YMR <- read.csv("/RNASEQ/YMR.matrix",sep="\t",
                stringsAsFactors = F,header = T)

YMR_1 <- as.matrix(YMR[,-1])

library(preprocessCore)
YMR_1<-normalize.quantiles(YMR_1)

YMR[,2:7] <- YMR_1[,1:6]

m.max = apply(YMR[,2:7],1,max)
quantile(m.var, probs = seq(0, 1, 0.01))
YMR = YMR[which(m.max>=1),]


write.table(HN1,"/RNASEQ/HN1_TPM_express.csv",
            sep = ",",col.names = T,row.names = F,quote = F)
write.table(YMR,"/RNASEQ/YMR_TPM_express.csv",
            sep = ",",col.names = T,row.names = F,quote = F)

```

### 1.2.2 Threshold 

```R
single_gene <- read.csv("/three_Papaver_species_protein/SingleCopyOrthogroups.txt",
                        sep="\t",header = F,stringsAsFactors = F)
gene_list <- read.csv("/three_Papaver_species_protein/Orthogroups.csv",
                      sep = "\t",stringsAsFactors = F)
gene_list <- gene_list[gene_list$X %in% single_gene[,1], ]

HN1 <- read.csv("/RNASEQ/HN1_TPM_express.csv",stringsAsFactors = F)
YMR <- read.csv("/RNASEQ/YMR_TPM_express.csv",stringsAsFactors = F)

HN1_YMR_expression <- gene_list[,c(3,4)]

HN1_YMR_expression <- merge(HN1_YMR_expression,HN1,by.x="HN1.pep",by.y="V2")
colnames(HN1_YMR_expression)[3:8] <- paste0("HN1_",colnames(HN1_YMR_expression)[3:8])

HN1_YMR_expression <- merge(HN1_YMR_expression,YMR,by.x="YMR.pep",by.y="V2")
colnames(HN1_YMR_expression)[9:14] <- paste0("YMR_",colnames(HN1_YMR_expression)[9:14])

HN1_YMR_expression[,15] <- NA

for (i in 1:dim(HN1_YMR_expression)[1]) {
  HN1_YMR_expression[i,15] <- cor(as.numeric(HN1_YMR_expression[i,3:8]),as.numeric(HN1_YMR_expression[i,9:14]),
                                  method = "spearman")
}  
mean(HN1_YMR_expression[,15])
##0.456217

```



### 1.2.3 Spearman correlation coefficients for duplicates

```R
retain_data <- read.csv("/HN1_YMR_alignment_final_rdup_chr.csv",
                        stringsAsFactors = F,header = F)
colnames(retain_data) <- c("YMR","subgenome1","subgenome2","chr")

both_retain <- retain_data[retain_data$subgenome1 !="",]
both_retain <- both_retain[both_retain$subgenome2 !="",]
both_retain <- both_retain[,-4]

both_retain[,4:6] <-NA
colnames(both_retain)[4:5] <- c("YMR_1","YMR_2")

HN1 <- read.csv("/RNASEQ/HN1_TPM_express.csv",stringsAsFactors = F)
YMR <- read.csv("/RNASEQ/YMR_TPM_express.csv",stringsAsFactors = F)


for (i in 1:dim(both_retain)[1]) {
  both_retain[i,4] <- cor(as.numeric(HN1[HN1$V2 == both_retain[i,2],][2:7]),
                          as.numeric(YMR[YMR$V2 == both_retain[i,1],][2:7]),
                          method = "spearman")
  
  both_retain[i,5] <- cor(as.numeric(HN1[HN1$V2 == both_retain[i,3],][2:7]),
                          as.numeric(YMR[YMR$V2 == both_retain[i,1],][2:7]),
                          method = "spearman")
  
}

both_retain <- na.omit(both_retain)

write.table(both_retain,"both_retain_cor.csv",
            sep = ",",col.names = T,row.names = F,quote = F)

```



### 1.2.4 Classify

```R
both_retain <- read.csv("both_retain_cor.csv",stringsAsFactors = F)
both_retain[,7] <- NA
colnames(both_retain)[7] <- "Type"

threshold_con <- 0.456217
threshold_spe <- 0.456217
for (i in 1:dim(both_retain)[1]) {
  if(both_retain[i,4] >= threshold_con && both_retain[i,5] >= threshold_con){
    both_retain[i,7] <- "Both_Conservation"
  }else if(both_retain[i,4] < threshold_spe && both_retain[i,5] >= threshold_con){
    both_retain[i,7] <- "Specialization_1"
  }else if(both_retain[i,4] >= threshold_con && both_retain[i,5] < threshold_spe){
    both_retain[i,7] <- "Specialization_2"
  }else if(both_retain[i,4] < threshold_spe && both_retain[i,5] < threshold_spe){
    both_retain[i,7] <- "Both_Specialization"
  }else{
    both_retain[i,7] <- "others"
  }
}

result_data <- as.data.frame(table(both_retain[,7]))

write.table(both_retain,paste0("/COR/","all",".csv"),
            sep = ",",col.names = T,row.names = F,quote = F)

for (i in result_data[,1]) {
  write.table(both_retain[both_retain$Type == i,],paste0("/COR/",i,".csv"),
              sep = ",",col.names = T,row.names = F,quote = F)
}


```



## 1.3 Singletons

### 1.3.1 Spearman correlation coefficients for singletons

```R
retain_data <- read.csv("/HN1_YMR_alignment_final_rdup_chr.csv",
                        stringsAsFactors = F,header = F)
retain_data <- retain_data[,-4]
colnames(retain_data) <- c("YMR","subgenome1","subgenome2")

both_retain <- retain_data[retain_data$subgenome1 !="",]
both_retain <- both_retain[both_retain$subgenome2 !="",]


subgenome1_singleton <- retain_data[retain_data$subgenome1 !="",]
subgenome1_singleton <- subgenome1_singleton[!subgenome1_singleton[,1] %in% both_retain[,1],]


subgenome2_singleton <- retain_data[retain_data$subgenome2 !="",]
subgenome2_singleton <- subgenome2_singleton[!subgenome2_singleton[,1] %in% both_retain[,1],]


subgenome1_singleton[,4] <-NA
colnames(subgenome1_singleton)[4] <- c("cor")

subgenome2_singleton[,4] <-NA
colnames(subgenome2_singleton)[4] <- c("cor")

HN1 <- read.csv("/RNASEQ/HN1_TPM_express.csv",stringsAsFactors = F)
YMR <- read.csv("/RNASEQ/YMR_TPM_express.csv",stringsAsFactors = F)


for (i in 1:dim(subgenome1_singleton)[1]) {
  subgenome1_singleton[i,4] <- cor(as.numeric(HN1[HN1$V2 == subgenome1_singleton[i,2],][2:7]),
                                   as.numeric(YMR[YMR$V2 == subgenome1_singleton[i,1],][2:7]),
                                   method = "spearman")
  
}

for (i in 1:dim(subgenome2_singleton)[1]) {
  subgenome2_singleton[i,4] <- cor(as.numeric(HN1[HN1$V2 == subgenome2_singleton[i,3],][2:7]),
                                   as.numeric(YMR[YMR$V2 == subgenome2_singleton[i,1],][2:7]),
                                   method = "spearman")
  
}

subgenome1_singleton <- na.omit(subgenome1_singleton)

subgenome2_singleton <- na.omit(subgenome2_singleton)

write.table(subgenome1_singleton,"/subgenome1_singleton_cor.csv",
            sep = ",",col.names = T,row.names = F,quote = F)
write.table(subgenome2_singleton,"/subgenome2_singleton_cor.csv",
            sep = ",",col.names = T,row.names = F,quote = F)

```



### 1.3.2 Classify

```R
subgenome1 <- read.csv("/subgenome1_singleton_cor.csv",stringsAsFactors = F)
subgenome1[,5] <- NA
colnames(subgenome1)[5] <- "Type"

threshold_con <- 0.456217
threshold_spe <- 0.456217
for (i in 1:dim(subgenome1)[1]) {
  if(subgenome1[i,4] >= threshold_con){
    subgenome1[i,5] <- "Conservation"
  }else if(subgenome1[i,4] < threshold_spe){
    subgenome1[i,5] <- "Specialization"
  }
}

result_data <- as.data.frame(table(subgenome1[,5]))

write.table(subgenome1,paste0("/COR/singleton/","subgenome1_all",".csv"),
            sep = ",",col.names = T,row.names = F,quote = F)

for (i in result_data[,1]) {
  write.table(subgenome1[subgenome1$Type == i,],paste0("/COR/singleton/subgenome1_",i,".csv"),
              sep = ",",col.names = T,row.names = F,quote = F)
}


subgenome2 <- read.csv("/subgenome2_singleton_cor.csv",stringsAsFactors = F)
subgenome2[,5] <- NA
colnames(subgenome2)[5] <- "Type"

threshold_con <- 0.456217
threshold_spe <- 0.456217
for (i in 1:dim(subgenome2)[1]) {
  if(subgenome2[i,4] >= threshold_con){
    subgenome2[i,5] <- "Conservation"
  }else if(subgenome2[i,4] < threshold_spe){
    subgenome2[i,5] <- "Specialization"
  }
}

result_data <- as.data.frame(table(subgenome2[,5]))

write.table(subgenome2,paste0("/COR/singleton/","subgenome2_all",".csv"),
            sep = ",",col.names = T,row.names = F,quote = F)

for (i in result_data[,1]) {
  write.table(subgenome2[subgenome2$Type == i,],paste0("/COR/singleton/subgenome2_",i,".csv"),
              sep = ",",col.names = T,row.names = F,quote = F)
}

```



## 1.4 Identify tissue-specific genes

```R
YMR <- read.csv("RNASEQ/YMR.matrix",sep="\t",
                stringsAsFactors = F,header = T)
YMR <- YMR[,-1]
YMR <- YMR[,c(7,1:6)]
YMR_1 <- as.matrix(YMR[,-1])

library(preprocessCore)
YMR_1<-normalize.quantiles(YMR_1)

YMR[,2:7] <- YMR_1[,1:6]
library("RNentropy")
data_matrix <- diag(x = 1, 6, 6)
Results <- RN_calc(YMR[1:dim(YMR)[1],], data_matrix)
Results <- RN_select(Results)
YMR <- YMR[YMR[,1] %in% Results$selected[,1],]


tissue_specific <- function(expression_vec){
  for (i in 1:6) {
    if(expression_vec[i] > 1 & expression_vec[i] > max(expression_vec[-i]) ){
      if(min(scale(expression_vec)) > -1.2){
        return(i)}
      
    }
  }
  return(0)
}

YMR[,8] <-  apply(YMR[,2:7],1,tissue_specific)
YMR_tissue_specific <- YMR[YMR[,8] != 0,]

for (i in 1:dim(YMR_tissue_specific)[1]) {
  YMR_tissue_specific[i,YMR_tissue_specific[i,8]+1] <- 1
  YMR_tissue_specific[i,-c(1,YMR_tissue_specific[i,8]+1,8)] <- 0
}

HN1 <- read.csv("/RNASEQ/HN1.matrix",sep="\t",
                stringsAsFactors = F,header = T)

HN1_1 <- as.matrix(HN1[,-1])

library(preprocessCore)
HN1_1<-normalize.quantiles(HN1_1)

HN1[,2:7] <- HN1_1[,1:6]


data_matrix <- diag(x = 1, 6, 6)
Results <- RN_calc(HN1[1:dim(HN1)[1],], data_matrix)
Results <- RN_select(Results)
HN1 <- HN1[HN1[,1] %in% Results$selected[,1],]

HN1[,8] <-  apply(HN1[,2:7],1,tissue_specific)
HN1_tissue_specific <- HN1[HN1[,8] != 0,]

for (i in 1:dim(HN1_tissue_specific)[1]) {
  HN1_tissue_specific[i,HN1_tissue_specific[i,8]+1] <- 1
  HN1_tissue_specific[i,-c(1,HN1_tissue_specific[i,8]+1,8)] <- 0
}


YMR <- YMR[YMR[,8] != 0,]
YMR <- YMR[YMR[,1] %in% YMR_tissue_specific[,1],]
HN1 <- HN1[HN1[,8] != 0,]
HN1 <- HN1[HN1[,1] %in% HN1_tissue_specific[,1],]

YMR <- YMR[order(YMR$V8),]
HN1 <- HN1[order(HN1$V8),]

write.csv(HN1,"HN1_tissue_specific.csv")
write.csv(YMR,"YMR_tissue_specific.csv")
library(pheatmap)

pheatmap(as.matrix(YMR[order(YMR$V8),][,c(-1,-8)]),show_rownames =FALSE,scale="row",
         cluster_rows = FALSE,cluster_cols = FALSE,main="YMR")
pheatmap(as.matrix(HN1[order(HN1$V8),][,c(-1,-8)]),show_rownames =FALSE,scale="row",
         cluster_rows = FALSE,cluster_cols = FALSE,main="HN1")


```

