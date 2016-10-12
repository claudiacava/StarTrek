#' @title Get human KEGG pathway data.
#' @description getKEGGdata creates a data frame with human KEGG pathway. Columns are the pathways and rows the genes inside those pathway 
#' @param KEGG_path  variable
#' @export
#' @importFrom KEGGREST keggList keggGet
#' @importFrom org.Hs.eg.db org.Hs.egSYMBOL2EG
#' @importFrom AnnotationDbi mappedkeys as.list
#' @return dataframe with human pathway data
#' @examples
#' path<-getKEGGdata(KEGG_path="Ener_met")
getKEGGdata<-function(KEGG_path){
if (KEGG_path=="Carb_met") {
  c<-proc_path(pathwayKEGG)
  a<-c[[2]]
}
  if (KEGG_path=="Ener_met") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Lip_met") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Amn_met") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Gly_bio_met") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Cof_vit_met") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Transcript") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Transl") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Fold_degr") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Repl_repair") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="sign_transd") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="sign_mol_int") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="Transp_cat") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="cell_grow_d") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="cell_comm") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="imm_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="end_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }
  if (KEGG_path=="circ_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  } 
  if (KEGG_path=="dig_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  } 
  if (KEGG_path=="exc_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  }  
  if (KEGG_path=="nerv_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  } 
  if (KEGG_path=="sens_syst") {
    c<-proc_path(pathwayKEGG)
    a<-c[[2]]
  } 
if (KEGG_path=="KEGG_path") {
  pathways.list <- keggList("pathway", "hsa")## returns the list of human pathways
pathway.codes <- sub("path:", "", names(pathways.list))
pathways.list<-list(pathways.list)
pathways.list<-pathways.list[lapply(pathways.list,length)!=0] 
a<-do.call("cbind", pathways.list)
}
pathway.codes<-c[[1]]
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
#' organism="Saccharomyces_cerevisiae"
#' netw<-getNETdata(network="SHpd",organism)
getNETdata<-function(network,organism=NULL){
  org_shar_pro<-SpidermiRquery_species(species)
  if (is.null(organism)) {
  net_shar_prot<-SpidermiRquery_spec_networks(organismID = org_shar_pro[6,],network)
  out_net_shar_pro<-SpidermiRdownload_net(net_shar_prot)
  geneSymb_net_shar_pro<-SpidermiRprepare_NET(organismID = org_shar_pro[6,],data = out_net_shar_pro)
  }
  if( !is.null(organism) ){
    net_shar_prot<-SpidermiRquery_spec_networks(organismID = org_shar_pro[9,],network)
    out_net_shar_pro<-SpidermiRdownload_net(net_shar_prot)
    geneSymb_net_shar_pro<-SpidermiRprepare_NET(organismID = org_shar_pro[9,],data = out_net_shar_pro)
}
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





