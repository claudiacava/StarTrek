#' @title Get human KEGG pathway data.
#' @description getKEGGdata creates a data frame with human KEGG pathway. Columns are the pathways and rows the genes inside those pathway 
#' @param KEGG_path  variable
#' @export
#' @importFrom KEGGREST keggList keggGet
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL2EG
#' @importFrom AnnotationDbi mappedkeys as.list
#' @return dataframe with human pathway data
#' @examples
#' path<-getKEGGdata(KEGG_path="Carb_met")
getKEGGdata<-function(KEGG_path=NULL){

if (KEGG_path=="Carb_met") {
  pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
mer<-select_path_carb(Carbohydrate)
common<-intersect(pathways.list,mer)
lo<-list()
for (i in 1:length(pathways.list)){
if (length(intersect(pathways.list[[i]],common)!=0)){
lo[[i]]<-pathways.list[[i]]
names(lo)[[i]]<-names(pathways.list)[[i]]
}
}
pathways.list<-lo[lapply(lo,length)!=0] 
pathway.codes <- sub("path:", "", names(pathways.list))
a<-do.call("rbind", pathways.list)
}

  if (KEGG_path=="Ener_met") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_en(Energy)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }

  
  if (KEGG_path=="Lip_met") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_lip(Lipid)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  
  if (KEGG_path=="Amn_met") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_amn(Aminoacid)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  
  if (KEGG_path=="Gly_bio_met") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_gly(Glybio_met)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="Cof_vit_met") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_cofa(Cofa_vita_met)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="Transcript") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_transc(Transcription)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  
  
  if (KEGG_path=="Transl") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_transl(Translation)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  
  if (KEGG_path=="Fold_degr") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_fold(Folding_sorting_and_degradation)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="Repl_repair") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_repl(Replication_and_repair)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="sign_transd") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_sign(Signal_transduction)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="sign_mol_int") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_sign_mol(Signaling_molecules_and_interaction)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="Transp_cat") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_transp_ca(Transport_and_catabolism)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="cell_grow_d") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_cell_grow(Cell_growth_and_death)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="cell_comm") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_cell_comm(Cellular_community)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="imm_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_imm_syst(Immune_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="end_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_end_syst(Endocrine_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }
  
  if (KEGG_path=="circ_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_circ_syst(Circulatory_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  } 
  
  if (KEGG_path=="dig_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_dig_syst(Digestive_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  } 
  
  if (KEGG_path=="exc_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_exc_syst(Excretory_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  }  
  
  
  if (KEGG_path=="nerv_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_ner_syst(Nervous_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  } 
  
  if (KEGG_path=="sens_syst") {
    pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
    mer<-select_path_sens_syst(Sensory_system)
    common<-intersect(pathways.list,mer)
    lo<-list()
    for (i in 1:length(pathways.list)){
      if (length(intersect(pathways.list[[i]],common)!=0)){
        lo[[i]]<-pathways.list[[i]]
        names(lo)[[i]]<-names(pathways.list)[[i]]
      }
    }
    pathways.list<-lo[lapply(lo,length)!=0] 
    pathway.codes <- sub("path:", "", names(pathways.list))
    a<-do.call("rbind", pathways.list)
  } 
  
if (KEGG_path=="KEGG_path") {
  pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
pathway.codes <- sub("path:", "", names(pathways.list))
pathways.list<-list(pathways.list)
pathways.list<-pathways.list[lapply(pathways.list,length)!=0] 
a<-do.call("cbind", pathways.list)
}

genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             pw[[1]]$GENE[c(TRUE, FALSE)]
                           })
x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
top3 <- matrix(0, length(xx), length(genes.by.pathway))
rownames(top3) <- names(xx)
colnames(top3)<- names(genes.by.pathway)
for (j in  1:length(xx)){
  for (k in  1:length(genes.by.pathway)){
    if (length(intersect(xx[[j]],genes.by.pathway[[k]])!=0)){
      
      top3[j,k]<-names(xx[j]) 
    }
  }
}
top3[top3 == 0] <- " "



#a<-data.frame(pathways.list)
#i <- sapply(a, is.factor)
#a[i] <- lapply(a[i], as.character)
rownames(a)<-sub("path:","",rownames(a))
PROVA<-top3
ti=list()
for( i in 1:ncol(PROVA)) {
  if (colnames(PROVA)[i]==rownames(a)[i]){
    colnames(PROVA)[i]<-a[i]
}
}
return(PROVA)
}


#' @title Get network data.
#' @description getNETdata creates a data frame with network data. 
#' Network category can be filtered among: physical interactions, co-localization, genetic interactions and shared protein domain.
#' @param network  variable. The user can use the following parameters 
#' based on the network types to be used. PHint for Physical_interactions,
#' COloc for Co-localization, GENint for Genetic_interactions and
#' SHpd for Shared_protein_domains
#' @export
#' @importFrom SpidermiR SpidermiRquery_species SpidermiRquery_spec_networks SpidermiRdownload_net SpidermiRprepare_NET
#' @return dataframe with gene-gene (or protein-protein interactions)
#' @examples
#' netw<-getNETdata(network="SHpd")
getNETdata<-function(network){
  org_shar_pro<-SpidermiRquery_species(species)
  net_shar_prot<-SpidermiRquery_spec_networks(organismID = org_shar_pro[6,],network)
  out_net_shar_pro<-SpidermiRdownload_net(net_shar_prot)
  geneSymb_net_shar_pro<-SpidermiRprepare_NET(organismID = org_shar_pro[6,],data = out_net_shar_pro)
  ds_shar_pro<-do.call("rbind", geneSymb_net_shar_pro)
  data_shar_pro<-as.data.frame(ds_shar_pro[!duplicated(ds_shar_pro), ]) 
  sdc_shar_pro<-unlist(data_shar_pro$gene_symbolA,data_shar_pro$gene_symbolB)
  m_shar_pro<-c(data_shar_pro$gene_symbolA)
  m2_shar_pro<-c(data_shar_pro$gene_symbolB)
  ss_shar_pro<-cbind(m_shar_pro,m2_shar_pro)
  data_pr_shar_pro<-as.data.frame(ss_shar_pro[!duplicated(ss_shar_pro), ]) 
  colnames(data_pr_shar_pro) <- c("m_shar_pro", "m2_shar_pro")
return(data_pr_shar_pro)
}


#' @title Get TCGA data.
#' @description get_TCGAdata creates a data frame with network data. 
#' Network category can be filtered among: physical interactions, co-localization, genetic interactions and shared protein domain.
#' @param cancer cancer type, See TCGAbiolinks package
#' @param PlatformCancer platform, type See TCGAbiolinks package
#' @param tumour barcode samples with label for example tumour
#' @param normal barcode samples with label for example normal
#' @param patha directory name
#' @export
#' @importFrom SpidermiR SpidermiRquery_species SpidermiRquery_spec_networks SpidermiRdownload_net SpidermiRprepare_NET
#' @importFrom TCGAbiolinks GDCquery GDCdownload GDCprepare TCGAanalyze_Filtering
#' @return dataframe with gene-gene (or protein-protein interactions)
#' @examples
#' cancer <- "TCGA-BRCA"
#' PlatformCancer <- "Illumina HiSeq"
#' tumour<-c("TCGA-BH-A0DL-01A-11R-A115-07","TCGA-AO-A03P-01A-11R-A00Z-07")
#' normal<-c("TCGA-BH-A209-11A-42R-A157-07","TCGA-E9-A1N4-11A-33R-A14M-07") 
#' TCGA_matrix<-get_TCGAdata(cancer,PlatformCancer,tumour,normal,
#'                                        patha = "exampleData")
get_TCGAdata<-function(cancer,PlatformCancer,tumour,normal,patha = "exampleData"){
dataType <- "normalized_results"  
query <- GDCquery(project = cancer,
                  data.category = "Gene expression",
                  data.type = "Gene expression quantification",
                  platform = PlatformCancer, 
                  file.type  = dataType, 
                  barcode = c(tumour,normal),
                  legacy = TRUE)
GDCdownload(query,directory = patha)

dataAssy <- GDCprepare(query,  directory = patha, summarizedExperiment = FALSE)

dataFilt <- TCGAanalyze_Filtering(tabDF = dataAssy,
                                  method = "quantile", 
                                  qnt.cut =  0.25)  

colnames(dataFilt) <- gsub("normalized_count_","",colnames(dataFilt))







return(dataFilt)
}


