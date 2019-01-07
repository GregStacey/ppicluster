
if (!require(ontologyIndex)) {
  install.packages("ontologyIndex", repos='http://cran.us.r-project.org')
}
library(ontologyIndex)

if (!require(GOstats)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("GOstats")
}
library(GOstats)

if (!require("tidyverse")) {
  install.packages("tidyverse", repos='http://cran.us.r-project.org')
}
library(tidyverse)

if (!require("magrittr")) {
  install.packages("magrittr", repos='http://cran.us.r-project.org')
}
library(magrittr)

if (!require("flavin")) {
  install.packages("devtools")
  require(devtools)
  install_github("skinnider/flavin")
}
library(flavin)

if (!require("org.Hs.eg.db")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db")
}
library(org.Hs.eg.db)



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
slim = get_ontology("../data/goslim_generic.obo")
goa.slim = read_gpa("../data/goa_human.gpa",
                    filter.NOT = T, filter.evidence = "IPI",
                    ontology = slim, propagate = T)
goa.slim = goa.slim[goa.slim$UNIPROT %in% unqprots,]


# make a uniprot <--> ensembl map
fn.map = "../data/HUMAN_9606_idmapping.dat.gz"
map = read_tsv(fn.map,
               col_names = c("uniprot", "db", "id"))
ens = filter(map, db == "GeneID") %>% dplyr::select(-db)
ent = filter(map, db == 'UniProtKB-ID') %>% dplyr::select(-db)
ens_ent = left_join(ens, ent, by = 'uniprot') %>%
  #dplyr::select(-uniprot) %>%
  drop_na() %>%
  set_colnames(c("uniprot", "entrez", "gene"))


# read gene2go file in case GOstats fails
fn = "../data/gene2go"
gene2go = read_tsv(fn)
gene2go = gene2go[!gene2go$Evidence %in% "IEA", ] # filter by evidence
# filter to go slim
gene2go = gene2go[gene2go$GO_ID %in% goa.slim$GO.ID, ]


clusters$n.BP.sig = -1
# cluster set
for (ii in 1:length(allsets)) {
  I.clusters = which(clusters$set %in% allsets[ii])
  this.background = unique(unlist(strsplit(clusters$cluster[I.clusters], ";")))
  this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
  
  # complex
  for (jj in 1:length(I.clusters)) {
    print(paste("enrichment analysis ", allsets[ii], " complex", jj))
    
    this.complex = clusters$cluster[I.clusters[jj]]
    this.ids = unlist(strsplit(clusters$cluster[I.clusters[jj]], ";"))
    this.ids.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
    
    clusters$n.BP.sig[I.clusters[jj]] = NA
    if (length(this.ids.entrez)<2) {
      print("    cluster size <= 2...")
      next
    }
    
    # need this nonsense trycatch() because of an sqlite error
    # on computecanada
    #sumTable = NA
    #sumTable <- tryCatch(
    #  {
    #    param <- new("GOHyperGParams", geneIds=this.ids.entrez,
    #                 universeGeneIds=this.background.entrez,
    #                 annotation="org.Hs.eg.db", ontology="BP",pvalueCutoff=0.05,
    #                 conditional=FALSE, testDirection="over")
    #    hyp <- hyperGTest(param)
    #    sumTable <- summary(hyp)
    #  },
    #  error=function(cond) {

    # get all GO terms associated with background
    this.gene2go = gene2go[gene2go$GeneID %in% this.background.entrez,]
    
    # remove redundant entries
    this.gene2go = distinct(this.gene2go[,c("GeneID", "GO_ID")])
    this.goterms = unique(this.gene2go$GO_ID)
    nn = length(this.goterms)
    
    print(paste("    testing", nn, "go bp terms..."))
    
    sumTable = data.frame(GOBPID = character(nn), Pvalue = numeric(nn),
                          stringsAsFactors = F)
    for (kk in 1:nn) {
      go = this.goterms[kk]
      sumTable$GOBPID[kk] = go
      #print(sumTable$GOBPID[kk])
      
      # universe
      I.clusters = which(clusters$set %in% allsets[ii])
      this.background = unique(unlist(strsplit(clusters$cluster[I.clusters], ";")))
      this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
      n.universe = length(this.background.entrez)
      
      # universe hits
      n.u.hits = sum(this.gene2go$GeneID %in% this.background.entrez & this.gene2go$GO_ID %in% go)
      
      # set
      this.complex = clusters$cluster[I.clusters[jj]]
      this.ids = unlist(strsplit(clusters$cluster[I.clusters[jj]], ";"))
      this.ids.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
      n.set = length(this.ids.entrez)
      
      # set hits
      n.s.hits = sum(this.gene2go$GeneID %in% this.ids.entrez & this.gene2go$GO_ID %in% go)
      
      sumTable$Pvalue[kk] = phyper(n.s.hits-1, n.u.hits, n.universe-n.u.hits, n.set, lower.tail=FALSE)
    }
    sumTable = sumTable[sumTable$Pvalue<=.01, ]
    #},
    #warning=function(cond) {
    #  print("error")
    #},
    #finally={
    #  print("umm, finally")
    #}
    #)    
    if (is.data.frame(sumTable)) {
      clusters$n.BP.sig[I.clusters[jj]] = nrow(sumTable)
      print(paste("    number of sig terms = ", nrow(sumTable)))
    }
  }
  
  # write after each set in case of crash
  write_tsv(clusters, path = "../data/BP_clusters.txt")
}

# write
write_tsv(clusters, path = "../data/BP_clusters.txt")



## analyze BP_clusters table
fn = "../data/BP_clusters_graham.txt"
bp = read_tsv(fn)
bp = bp[!bp$noise_mag == -1,]

# fraction of clusters with no support
nn = length(allsets)
df = data.frame(data_type=character(nn),
                noise_type=character(nn),
                noise_mag=numeric(nn),
                algorithm = character(nn),
                mean = numeric(nn),
                fraction0 = numeric(nn),
                nn = numeric(nn), stringsAsFactors = F)
for (ii in 1:length(allsets)) {
  tmp = unlist(strsplit(allsets[ii], ";"))
  I = bp$data_type%in%tmp[1] & bp$noise_type%in%tmp[2] & 
    bp$noise_mag==as.numeric(tmp[3]) & bp$algorithm%in%tmp[4] &
    bp$n.BP.sig>-1 & bp$noise_type %in% c("chrom","network_remove","network_add", "network_shuffle")
  
  df$data_type[ii] = tmp[1]
  df$noise_type[ii] = tmp[2]
  df$noise_mag[ii] = as.numeric(tmp[3])
  df$algorithm[ii] = tmp[4]
  df$mean[ii] = mean(bp$n.BP.sig[I], na.rm=T)
  df$fraction0[ii] = mean(bp$n.BP.sig[I]==0, na.rm=T)
  df$nn[ii] = sum(!is.na(bp$n.BP.sig[I]))
}
df = df[!df$noise_type %in% "-1",]


I = bp$noise_mag>-1 & bp$n.BP.sig>-1 & !is.na(log10(bp$n.BP.sig+1)) & bp$noise_mag<=1
# number of sig terms
ggplot(bp[I,], aes(jitter(x=noise_mag), y=log10(n.BP.sig+1))) + 
  geom_point(alpha=0.05) + geom_smooth() + facet_grid(noise_type~algorithm)
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/GOvalidation_avgNterms.png")

# fraction 0
ggplot(df, aes(x=jitter(noise_mag), y=fraction0)) + 
  geom_point(alpha=0.3) + facet_grid(algorithm~noise_type, scales="free") + 
  geom_smooth()
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/GOvalidation_frac0terms.png")


