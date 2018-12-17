
if (!library(ontologyIndex)) {
  install.packages("ontologyIndex")
}
library(ontologyIndex)

if (!library(GOstats)) {
  install.packages("GOstats")
}
library(GOstats)

if (!library("ggplot2")) {
  install.packages("ggplot2")
}
require(ggplot2)

if (!library("tidyverse")) {
  install.packages("tidyverse")
}
library(tidyverse)

if (!library("magrittr")) {
  install.packages("magrittr")
}
library(magrittr)

if (!library("flavin")) {
  install.packages("flavin")
}
library(flavin)



fnsave = "../data/clusters_coexp.Rda"
load(fnsave)

# # load clusters
# fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters.txt"
# clusters = read.csv2(fn, sep="\t", quote="", stringsAsFactors = F)
# clusters = clusters[,!names(clusters) %in% "X"]
# clusters$noise_mag = as.numeric(clusters$noise_mag)
# 
# # add real CORUM complexes as positive control
# fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/allComplexes.txt"
# corum = read.table(fn, sep="\t", stringsAsFactors = F, quote="", header = T)
# corum2add = data.frame(data_type=character(nrow(corum)),
#                        noise_type=character(nrow(corum)),
#                        noise_mag=numeric(nrow(corum)),
#                        algorithm=character(nrow(corum)),
#                        cluster=corum$subunits.UniProt.IDs.,
#                        stringsAsFactors = F)
# corum2add$data_type = "corum_real"
# corum2add$noise_type = "-1"
# corum2add$noise_mag = "-1"
# corum2add$algorithm = "-1"
# clusters = rbind(clusters, corum2add)

# unique IDs
unqprots = unique(unlist(strsplit(clusters$cluster, ";")))

# assign each cluster to a set (data_type, noise_type, noise_mag, algorithm)
cols = c("data_type", "noise_type", "noise_mag", "algorithm")
clusters$set = do.call(paste, c(clusters[cols], sep=";"))
allsets = unique(clusters$set)


# read GO annotations
# ontology = get_ontology("/Users/Mercy/Academics/Foster/PCP-SILAC/Validation_R_skinnider_2/data/db/go-basic.obo")
# goa = read_gpa("/Users/Mercy/Academics/Foster/Manuscripts/GoldStandards/Data/goa_human.gpa",
#                filter.NOT = T, filter.evidence = "IPI",
#                ontology = ontology, propagate = T)
# goa = goa[goa$UNIPROT %in% unqprots,]
# read go slim
slim = get_ontology("../data/goslim_pir.obo")
goa.slim = read_gpa("/Users/Mercy/Academics/Foster/Manuscripts/GoldStandards/Data/goa_human.gpa",
               filter.NOT = T, filter.evidence = "IPI",
               ontology = slim, propagate = T)
goa.slim = goa.slim[goa.slim$UNIPROT %in% unqprots,]


# make a uniprot <--> ensembl map
fn.map = "/Users/Mercy/Academics/Foster/Data/dbs/HUMAN_9606_idmapping.dat.gz"
map = read_tsv(fn.map,
               col_names = c("uniprot", "db", "id"))
ens = filter(map, db == "GeneID") %>% dplyr::select(-db)
ent = filter(map, db == 'UniProtKB-ID') %>% dplyr::select(-db)
ens_ent = left_join(ens, ent, by = 'uniprot') %>%
  #dplyr::select(-uniprot) %>%
  drop_na() %>%
  set_colnames(c("uniprot", "entrez", "gene"))


clusters$n.sig.annot = numeric(nrow(clusters))
# cluster set
for (ii in 1:length(allsets)) {
  I.clusters = which(clusters$set %in% allsets[ii])
  this.background = unique(unlist(strsplit(clusters$cluster[I.clusters], ";")))
  this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
  
  # complex
  for (jj in 1:length(I.clusters)) {
    this.complex = clusters$cluster[I.clusters[jj]]
    this.ids = unlist(strsplit(clusters$cluster[I.clusters[jj]], ";"))
    this.ids.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
    
    clusters$n.sig.annot[I.clusters[jj]] = NA
    if (length(this.ids.entrez)<2) next
    
    param <- new("GOHyperGParams", geneIds=this.ids.entrez,
                 universeGeneIds=this.background.entrez,
                 annotation="org.Hs.eg.db", ontology="BP",pvalueCutoff=0.05,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp, categorySize=30)
    
    # filter to go slim
    sumTable = sumTable[sumTable$GOBPID %in% goa.slim$GO.ID,]
    
    clusters$n.sig.annot[I.clusters[jj]] = nrow(sumTable)
    print(paste("number of sig annotation terms = ", nrow(sumTable)))
  }
}

# write
write_tsv(clusters, path = "../data/test_clusters.txt")



# fn = "/Users/Mercy/Downloads/gene2go"
# gene2go = read_tsv(fn)
# 
# 
# # do a test (GO:0006412)
# go = "GO:0006412"
# # universe
# ii = 546
# I.clusters = which(clusters$set %in% allsets[ii])
# this.background = unique(unlist(strsplit(clusters$cluster[I.clusters], ";")))
# this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
# n.universe = length(this.background.entrez)
# 
# # universe hits
# n.u.hits = sum(gene2go$GeneID %in% this.background.entrez & gene2go$GO_ID %in% go)
# 
# # set
# jj = 70
# this.complex = clusters$cluster[I.clusters[jj]]
# this.ids = unlist(strsplit(clusters$cluster[I.clusters[jj]], ";"))
# this.ids.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
# n.set = length(this.ids.entrez)
# 
# # set hits
# n.s.hits = sum(gene2go$GeneID %in% this.ids.entrez & gene2go$GO_ID %in% go)
# 
# phyper(n.s.hits, n.u.hits, n.universe-n.u.hits, n.set, lower.tail=FALSE)
# 
# 
