#' @title Get human KEGG pathway data and network data in order to define the common gene.
#' @description path_net creates a list of network data for each human pathway. The network data will be generated when interacting genes belong to that pathway.  
#' @param net_type  network data as provided by getNETdata
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a list of network data for each pathway (interacting genes belong to that pathway)
#' @examples
#' \dontrun{
#' lista_net<-path_net(pathway=path,net_type=netw)
#' }
path_net<-function(pathway,net_type){
    lista_int<-list()
    for (k in 1:ncol(pathway)){
      #k=1 
      print(paste(k,"PATHWAY",colnames(pathway)[k]))
      currentPathway_genes<-pathway[,k]
      common1 <- intersect( net_type$m_shar_pro, currentPathway_genes)
      common2 <- intersect( net_type$m2_shar_pro, currentPathway_genes)
      if (length(common1)==0 & length(common2)==0 ){
        mago2<-character(length = 0)
      }
      if (length(common1)!=0 | length(common2)!=0 ){
        b=list()
        for (i in 1:length(common1)){
          x<-common1[i]
          n<-overlap(net_type,x,currentPathway_genes)
          b[[i]]<-n
        }
        v<-do.call("rbind", b)
        c=list()
        for (i in 1:length(common2)){
          x<-common1[i]
          n<-overlap(net_type,x,currentPathway_genes)
          c[[i]]<-n
        }
        v2<-do.call("rbind", b)
        mago<-rbind(v,v2)
        mago2<-mago[!duplicated(mago), ]
      }
      
      if (length(mago2)!=0){
        lista_int[[k]]<-mago2
      }
      if (length(mago2)==0){
        lista_int[[k]]<-"0"} 
    }   
    return(lista_int)
  }
  
  
  





#' @title Get human KEGG pathway data and network data in order to define the common gene.
#' @description list_path_net creates a list of interacting genes for each human pathway.   
#' @param net_type  network data as provided by getNETdata
#' @param pathway  pathway data as provided by getKEGGdata
#' @export
#' @return a list of genes for each pathway (interacting genes belong to that pathway)
#' @examples
#' \dontrun{
#' path<-getKEGGdata(KEGG_path="Transcript")
#' netw<-getNETdata(network="SHpd")
#' list_path<-list_path_net(net_type=netw,pathway=path)}
list_path_net<-function(net_type,pathway){
  i <- sapply(net_type, is.factor) 
  net_type[i] <- lapply(net_type[i], as.character)
  m<-c(net_type$m_shar_pro)
  m2<-c(net_type$m2_shar_pro)
  s<-c(m,m2)
  fr<- unique(s)
  n<-as.data.frame(fr)
  i <- sapply(n, is.factor) 
  n[i] <- lapply(n[i], as.character)
  matrice <- matrix(0, nrow(n), ncol(pathway))
  rownames(matrice) <- n[,1]
  colnames(matrice)<- colnames(pathway)
  for (j in  1:nrow(n)){
   # k=1 
    for (k in  1:ncol(pathway)){
      if (length(intersect(n[j,],pathway[,k])!=0)){
        matrice[j,k]<-n[j,]
      }
    }
  }
  return(matrice)
}



#' @title For TCGA data get human pathway data and creates a matrix with the average of genes for each pathway.
#' @description average creates a matrix with a summarized value for each pathway  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' score_mean<-average(TCGA_matrix,path)
average<-function(dataFilt,pathway){
  DataMatrix<-dataFilt
  dataFilt[ , "new.col"] <- gsub("\\|.*", "", rownames(dataFilt))
  DataMatrix<-dataFilt[which(dataFilt$new.col!="?"),]
  DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
  rownames(DataMatrix)<-DataMatrix$new.col
  DataMatrix$new.col<-NULL
#DataMatrix <- cbind(Data_CANCER_normUQ_filt_normal, Data_CANCER_mRNA_tumor_HighStage)
#DataMatrix <- Data_CANCER_mRNA_tumor_HighStage
PEAmatrix <- matrix( 0, ncol(pathway),ncol(DataMatrix))
rownames(PEAmatrix) <- colnames(pathway)
colnames(PEAmatrix) <-  colnames(DataMatrix)
listIPA_pathways<-colnames(pathway)
for ( k in 1: nrow(PEAmatrix)){
  #k=1
  currentPathway <- colnames(pathway)[k]
  currentPathway_genes_list_common <- intersect(rownames(DataMatrix), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- DataMatrix[currentPathway_genes_list_common,]
  SumGenes <- colSums(currentPathway_genes_list_commonMatrix)
  AverageGenes <- SumGenes / length(currentPathway_genes_list_common)
  PEAmatrix[k,] <- AverageGenes
}
return(PEAmatrix)
}




#' @title For TCGA data get human pathway data and creates a matrix with the median of genes for each pathway.
#' @description score_median creates a matrix with a summarized value for each pathway  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' med<-score_median(TCGA_matrix,path)
score_median<-function(dataFilt,pathway){
  DataMatrix<-dataFilt
  DataMatrix[ , "new.col"] <- gsub("\\|.*", "", rownames(DataMatrix))
  DataMatrix<-DataMatrix[which(DataMatrix$new.col!="?"),]
  DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
  rownames(DataMatrix)<-DataMatrix$new.col
  DataMatrix$new.col<-NULL
  #DataMatrix <- cbind(Data_CANCER_normUQ_filt_normal, Data_CANCER_mRNA_tumor_HighStage)
  #DataMatrix <- Data_CANCER_mRNA_tumor_HighStage
  PEAmatrix <- matrix( 0, ncol(pathway),ncol(DataMatrix))
  rownames(PEAmatrix) <- colnames(pathway)
  colnames(PEAmatrix) <-  colnames(DataMatrix)
  listIPA_pathways<-colnames(pathway)
  for ( k in 1: nrow(PEAmatrix)){
    #k=1
    #currentPathway <- colnames(pathway)[k]
    currentPathway_genes_list_common <- intersect(rownames(DataMatrix), currentPathway_genes<-pathway[,k])
    currentPathway_genes_list_commonMatrix <- DataMatrix[currentPathway_genes_list_common,]
    mu<-list()
    for (i in 1:ncol(currentPathway_genes_list_commonMatrix)){
    med<-median(currentPathway_genes_list_commonMatrix[,i],na.rm = TRUE)
    mu[[i]]<-med
    }
    v<-do.call("cbind", mu)
    PEAmatrix[k,] <- v
  }
  return(PEAmatrix)
}


#' @title For TCGA data get human pathway data and creates a measure of cross-talk among pathways 
#' @description tau_dist creates a matrix with tau distance for pairwise pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @importFrom bioDist tau.dist
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' score_tau_dist<-tau_dist(TCGA_matrix,path)
tau_dist <- function(dataFilt,pathway){
  PEAmatrix<-average(dataFilt,pathway)
  datam = as.matrix(PEAmatrix)
  s1 = tau.dist(datam)	
  s2 = as.matrix(s1)
  return(s2)
}



#' @title For TCGA data get human pathway data and creates a measure of cross-talk among pathways 
#' @description euc_dist_crtlk creates a matrix with euclidean distance for pairwise pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' score_euc_dist<-euc_dist_crtlk(TCGA_matrix,path)
euc_dist_crtlk <- function(dataFilt,pathway){
  PEAmatrix<-average(dataFilt,pathway)
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  euc_dist<-cbind(df,ma_d) # inserisco label con le relazioni tra i pathway
  return(euc_dist)
}


#' @title For TCGA data get human pathway data and creates a measure of standard deviations among pathways 
#' @description st_dv creates a matrix with standard deviation for pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' stand_dev<-st_dv(dataFilt=TCGA_matrix,pathway=path)
st_dv<-function(dataFilt,pathway){
DataMatrix<-dataFilt
dataFilt[ , "new.col"] <- gsub("\\|.*", "", rownames(dataFilt))
DataMatrix<-dataFilt[which(dataFilt$new.col!="?"),]
DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
rownames(DataMatrix)<-DataMatrix$new.col
DataMatrix$new.col<-NULL
PEAmatrix_sd <- matrix( 0, ncol(pathway),ncol(DataMatrix))
rownames(PEAmatrix_sd) <- colnames(pathway)
colnames(PEAmatrix_sd) <-  colnames(DataMatrix)
for ( k in 1: nrow(PEAmatrix_sd)){
  currentPathway <- colnames(pathway)[k]
  currentPathway_genes_list_common <- intersect( rownames(DataMatrix), currentPathway_genes<-pathway[,k])
  currentPathway_genes_list_commonMatrix <- DataMatrix[currentPathway_genes_list_common,]
  stdev<-apply(currentPathway_genes_list_commonMatrix,2,sd) #deviazione standard dei pathway
  PEAmatrix_sd[k,] <- stdev}
return(PEAmatrix_sd)
}



#' @title For TCGA data get human pathway data and creates a measure of standard deviations among pathways 
#' @description st_dv creates a matrix with standard deviation for pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' cross_talk_st_dv<-dev_std_crtlk(dataFilt=TCGA_matrix,pathway=path)
dev_std_crtlk<-function(dataFilt,pathway){
PEAmatrix_sd<-st_dv(dataFilt,pathway)
df=combn(rownames(PEAmatrix_sd),2) 
df=t(df)

ma<-matrix(0,nrow(df),ncol(PEAmatrix_sd)) # creo matrix che conterr? le somme delle dev st
colnames(ma)<-colnames(PEAmatrix_sd) # colnames contiene il nome dei pazienti

for ( p in 1: ncol(PEAmatrix_sd)){ # per ogni paziente
  patients <- (PEAmatrix_sd)[,p] 
  out <- apply(df, 1, function(x) sum(patients[x])) # calcolo somma delle dev standard tra le possibili combinazioni
  ma[,p]<-out
}
somma_sd<-cbind(df,ma) 
return(somma_sd)
}


#' @title For TCGA data get human pathway data and creates a measure of discriminating score among pathways 
#' @description st_dv creates a matrix with standard deviation for pathways  
#' @param dataFilt TCGA matrix
#' @param pathway pathway data
#' @export
#' @return a matrix value for each pathway 
#' @examples
#' \dontrun{
#' path<-getKEGGdata(KEGG_path="Transcript")
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' cross_talk_st_dv<-ds_score_crtlk(dataFilt=TCGA_matrix,pathway=path)}
ds_score_crtlk<-function(dataFilt,pathway){
  PEAmatrix<-average(dataFilt,pathway)
  #step 5 distance
  # EUCLIDEA DISTANCE
  df=combn(rownames(PEAmatrix),2) # possibili relazioni tra i pathway
  df=t(df)
  ma_d<-matrix(0,nrow(df),ncol(PEAmatrix)) # creo matrix che conterr? le distanze
  colnames(ma_d)<-colnames(PEAmatrix) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix)){ # per ogni paziente
    patients <- (PEAmatrix)[,p] 
    distance<-dist(patients) # calcolo distanza EUCLIDEA tra le possibile combinazioni
    ma_d[,p]<-distance
  }
  PEAmatrix_sd<-st_dv(dataFilt,pathway)
  df=combn(rownames(PEAmatrix_sd),2) 
  df=t(df)
  ma<-matrix(0,nrow(df),ncol(PEAmatrix_sd)) # creo matrix che conterr? le somme delle dev st
  colnames(ma)<-colnames(PEAmatrix_sd) # colnames conterr? il nome dei pazienti
  for ( p in 1: ncol(PEAmatrix_sd)){ # per ogni paziente
    patients <- (PEAmatrix_sd)[,p] 
    out <- apply(df, 1, function(x) sum(patients[x])) # calcolo somma delle dev standard tra le possibili combinazioni
    ma[,p]<-out
  }
  score<-ma_d/ma # discriminating score M1-M2/S1+S2
  score<- cbind(df,score)  
return(score)
}

#' @title Hamming distance  
#' @description hamm_dist creates a matrix with Hamming distance values for pathways  
#' @param list_path pathway data
#' @export
#' @importFrom e1071  hamming.distance
#' @return a symmetric matrix value with a distance among pathways 
#' @examples
#' path<-getKEGGdata(KEGG_path="Transcript")
#' netw<-getNETdata(network="SHpd")
#' list_path<-list_path_net(net_type=netw,pathway=path)
#' hamm_distance<-hamm_dist(list_path)
hamm_dist<-function(list_pat){
alt<-t(list_pat)
d<-hamming.distance(as.matrix(alt))
return(d)}


#' @title SVM classification for each feature
#' @description svm class creates a list with auc value  
#' @param TCGA_matrix gene expression matrix
#' @param nfs nfs split data into a training  and test set
#' @param tumour barcode samples for a class
#' @param normal barcode samples for another class
#' @export
#' @importFrom e1071 tune svm 
#' @importFrom ROCR prediction performance 
#' @importFrom  grDevices rainbow
#' @return a symmetric matrix value with a distance among pathways 
#' @examples
#' \dontrun{
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumo<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' norm<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_data<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
#' nf <- 60
#' res_class<-svm_classification(TCGA_matrix=TCGA_data,nfs=nf,normal=norm,tumour=tumo)
#' }
 
svm_classification<-function(TCGA_matrix,tumour,normal,nfs){
  #library("e1071")
  #library(ROCR)
  dataFilt<-TCGA_matrix
  DataMatrix<-dataFilt
  dataFilt[ , "new.col"] <- gsub("\\|.*", "", rownames(dataFilt))
  DataMatrix<-dataFilt[which(dataFilt$new.col!="?"),]
  DataMatrix <- subset(DataMatrix, !duplicated(DataMatrix$new.col)) 
  rownames(DataMatrix)<-DataMatrix$new.col
  DataMatrix$new.col<-NULL
  
  tDataMatrix<-as.data.frame(t(DataMatrix))

  tDataMatrix$Target<-0
  tum<-intersect(rownames(tDataMatrix),tumour)
  nor<-intersect(rownames(tDataMatrix),normal)
  tDataMatrix$
    
  Dataset_g1<-tDataMatrix[nor,]
  Dataset_g3<- tDataMatrix[tum,]
    
  
#training=read.table('C:/Users/UserInLab05/Desktop/trai.txt',header = TRUE)
#testset=read.table('C:/Users/UserInLab05/Desktop/test.txt',header = TRUE)

  #Dataset_g1 <- Data_CANCER_mRNA_tumor_HighStage[Data_CANCER_mRNA_tumor_HighStage$Target == 0, ]
#Dataset_g3 <- Data_CANCER_mRNA_tumor_HighStage[Data_CANCER_mRNA_tumor_HighStage$Target == 1, ]
  
tab_g1_training <- sample(Dataset_g1$ID,round(nrow(Dataset_g1) / 100 * nfs ))
tab_g3_training <- sample(Dataset_g3$ID,round(nrow(Dataset_g3) / 100 * nfs ))
tab_g1_testing <- as.factor(setdiff(Dataset_g1$ID,tab_g1_training))
tab_g3_testing <- as.factor(setdiff(Dataset_g3$ID,tab_g3_training))

FR<-intersect(Dataset_g1$ID,tab_g1_training)
rownames(Dataset_g1)<-Dataset_g1$ID
G1<-Dataset_g1[FR,]
FR1<-intersect(Dataset_g3$ID,tab_g3_training)
rownames(Dataset_g3)<-Dataset_g3$ID
G3<-Dataset_g3[FR1,]
training<-rbind(G1,G3)

F<-intersect(Dataset_g1$ID,tab_g1_testing)
rownames(Dataset_g1)<-Dataset_g1$ID

G1_testing<-Dataset_g1[FALSE,]

F1<-intersect(Dataset_g3$ID,tab_g3_testing)
rownames(Dataset_g3)<-Dataset_g3$ID
G3_testing<-Dataset_g3[F1,]

testing<-rbind(G1_testing,G3_testing)

## split data into a training (2/3) and test set (1/3)

#t<-training[,c(1,3)]
training[,2]<-NULL
x <- subset(training, select=-Target)
y <- training$Target
testing[,2]<-NULL
z<-subset(testing, select=-Target)

#z<-as.data.frame(z)

#colnames(z)<-"X2"

zi<-testing$Target

auc.df<-list()
svm_model_after_tune_COMPL<-list()
for( k in 2: ncol(training)){
  #print(k)
  svm_tune <- tune(svm, train.x=x, train.y=y, 
                   kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)),cross=10)
  svm_model_after_tune <- svm(Target ~ ., data=training[,c(1,k)], kernel="radial", cost=svm_tune$best.parameters[1], gamma=svm_tune$best.parameters[2],cross=10,probability = TRUE)
  #summary(svm_model_after_tune)

  j=k-1
  z2=z[,j]
  z3<-as.data.frame(z2)
  #colnames(z3)<-as.character(paste("X",j,sep = ""))
  colnames(z3)<-colnames(z)[j]
  #classifiersMatrix <- c(classifiersMatrix,svm_model_after_tune)
  pred <- predict(svm_model_after_tune,z3,decision.values=TRUE,cross=10)

  #a<-table(pred,zi)
  svm.roc <- prediction(attributes(pred)$decision.values, zi)
  svm.auc <- performance(svm.roc, 'tpr', 'fpr')

  perf <- performance(svm.roc, "auc")
  auc<-perf@y.values[[1]]
  
  auc.df[[j]]<- auc
  svm_model_after_tune_COMPL[[j]]<-svm_model_after_tune
  
  palette <- as.matrix(rainbow(ncol(z)))
  #print(j)
  if (j > 1) {
    plot(svm.auc,col=palette[j], add=TRUE)
  }
  else {
    plot(svm.auc, col=palette[j])
  }
  
  legend('bottomright', colnames(z), 
         lty=1, col=palette, bty='n', cex=.90,pch = 20,ncol=1)
  
  
}
names(auc.df) <- colnames(z)
return(auc.df)
}

