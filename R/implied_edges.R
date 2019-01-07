# Make a bar graph showing GO enrichment for two sets of co-complex edges: 
# those in the interactome and those not.


# TO DO
#   3. Finish coding "write df2write"
#   4. Edit figures to reflect noise_type/noise_mag and add try catch OR remove figures for later


if (!require("tidyverse")) {
  install.packages("tidyverse", repos='http://cran.us.r-project.org')
}
library(tidyverse)
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

if (!require(grex)) {
  install.packages("grex", repos='http://cran.us.r-project.org')
}


print("initializing...")
# load clusters
fnsave = "../data/clusters_wshuffle_coexp.Rda"
load(fnsave)
clusters$noise_mag = as.numeric(clusters$noise_mag)

# Janky command line argument handling:
this.args = as.numeric(commandArgs(trailingOnly = T))
if (length(this.args)==0) { # assume not called from the command line, run with noise_mag==0
  # reduce to just un-noised complexes
  clusters = clusters[clusters$noise_mag==0, ]
  noise_type = "all"
  noise_mag = 0
  
  print("using all noise_types with noise_mag=0")
} else if(length(this.args)==01) { # assume called correctly from implied_edges.sh
  # map command line argument (integer 1-27) to noise_type and noise_mag
  ia = mod(this.args-1, 8) + 1 # noise_mag index
  ib = ceiling(this.args / 8) # noise_type index
  print(paste("mapped", this.args, "to ia=", ia, "ib=", ib))
  
  noise_type = c("chrom","network_remove", "network_shuffle")[ib]
  noise_mag = sort(unique(clusters$noise_mag[clusters$noise_type %in% noise_type]))[ia+1]
  
  I = clusters$noise_type %in% noise_type & clusters$noise_mag %in% noise_mag
  clusters = clusters[I,]
} else {
  stop("bad command line argument")
}


print("get corum edges...")
# get corum edges
fn = "../data/allComplexes.txt"
corum = read_tsv(fn)
corum.edges = character(10^6)
cc = 0
for (ii in 1:nrow(corum)) {
  prots = unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";"))
  for (jj in 1:length(prots)) {
    for (kk in 1:length(prots)) {
      if (jj>=kk) next
      cc = cc+1
      tmp = sort(c(prots[jj], prots[kk]))
      corum.edges[cc] = paste(tmp,collapse="-")
    }
  }
}
corum.edges = corum.edges[1:cc]

print("get interactome edges...")
# get interactome edges
fn = "../data/interactomes.txt"
interactomes = read_tsv(fn)
interactomes$edges = apply(interactomes[,1:2],1,paste,collapse="-")

print("get cluster edges...")
# get cluster edges
nn = 10^6
edges = character(nn)
cc = 0
for (ii in 1:nrow(clusters)) {
  prots = unlist(strsplit(clusters$cluster[ii], ";"))
  
  for (jj in 1:length(prots)) {
    for (kk in 1:length(prots)) {
      if (jj>=kk) next
      cc = cc+1
      tmp = sort(c(prots[jj], prots[kk]))
      edges[cc] = paste(tmp,collapse="-")
    }
  }
}
edges = edges[1:cc]

# which edges are known?
edges.in.corum = edges %in% corum.edges
edges.in.interactome = edges %in% interactomes$edges



print("do annotation analysis...")
# now do annotation analysis

# unique IDs
unqprots = unique(unlist(strsplit(edges, "-")))

# read GO annotations
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
ens_ent$entrez = as.numeric(ens_ent$entrez)

# read gene2go file in case GOstats fails
fn = "../data/gene2go"
gene2go = read_tsv(fn)
gene2go = gene2go[!gene2go$Evidence %in% "IEA", ] # filter by evidence
# filter to go slim
gene2go = gene2go[gene2go$GO_ID %in% goa.slim$GO.ID, ]

# find the edges where both proteins are in gene2go
prots.single = unlist(strsplit(edges,"-"))
protA = prots.single[seq(from=1, to=length(prots.single), by=2)]
protB = prots.single[seq(from=2, to=length(prots.single), by=2)]
entA = ens_ent$entrez[match(protA, ens_ent$uniprot)]
entB = ens_ent$entrez[match(protB, ens_ent$uniprot)]
both.entrez = !is.na(entA) & !is.na(entB)
Igood = which(entA %in% gene2go$GeneID & entB %in% gene2go$GeneID)
# filter gene2go to just our proteins
gene2go_reduced = gene2go[gene2go$GeneID%in%entA | gene2go$GeneID%in%entB,]

# get terms for each entrez id
unq.entrez = unique(c(entA,entB))
ent.go = data.frame(ent=character(length(unq.entrez)),
                   go.bp=character(length(unq.entrez)),
                   go.cc=character(length(unq.entrez)),
                   go.mf=character(length(unq.entrez)),
                   stringsAsFactors = F)
for (ii in 1:length(unq.entrez)) {
  Ibp = gene2go_reduced$GeneID %in% unq.entrez[ii] & gene2go_reduced$Category %in% "Process"
  Icc = gene2go_reduced$GeneID %in% unq.entrez[ii] & gene2go_reduced$Category %in% "Component"
  Imf = gene2go_reduced$GeneID %in% unq.entrez[ii] & gene2go_reduced$Category %in% "Function"
  ent.go$ent[ii] = unq.entrez[ii]
  ent.go$go.bp[ii] = paste(gene2go_reduced$GO_ID[Ibp], collapse=";")
  ent.go$go.cc[ii] = paste(gene2go_reduced$GO_ID[Icc], collapse=";")
  ent.go$go.mf[ii] = paste(gene2go_reduced$GO_ID[Imf], collapse=";")
}

# for each edge, count the number of shared terms
N.terms = data.frame(BP=numeric(length(edges)),
                     CC=numeric(length(edges)),
                     MF=numeric(length(edges)), stringsAsFactors = F)
for (ii in 1:length(Igood)) {
  if (mod(ii,500)==1) {
    print(ii)
  }

  # count intersecting terms
  Ia = ent.go$ent %in% entA[Igood[ii]]
  Ib = ent.go$ent %in% entB[Igood[ii]]
  for (jj in 1:3) {
    go1 = unlist(strsplit(ent.go[Ia,jj+1],";"))
    go2 = unlist(strsplit(ent.go[Ib,jj+1],";"))
    N.terms[Igood[ii],jj] = length(intersect(go1,go2))
  }
}




print("do co-expression analysis...")
# now do co-expression analysis

# read gtex data
fn = "../data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
gtex = read.table(fn, sep="\t", quote="", stringsAsFactors = F, skip = 3)
gtex$ensembl = sapply(strsplit(gtex$V1, split=".", fixed=T), "[", 1)

# map uniprot_edges
data(gtexv7)
df.idmap = grex(gtexv7)
ens_mappable  = df.idmap$ensembl_id[df.idmap$uniprot_id %in% unqprots]
gtex = gtex[gtex$ensembl %in% ens_mappable,]

prots.single = unlist(strsplit(edges,"-"))
protA = prots.single[seq(from=1, to=length(prots.single), by=2)]
protB = prots.single[seq(from=2, to=length(prots.single), by=2)]
ensA = df.idmap$ensembl_id[match(protA, df.idmap$uniprot_id)]
ensB = df.idmap$ensembl_id[match(protB, df.idmap$uniprot_id)]

# dumb correction: make gtex columns numeric
gtex[,3:(ncol(gtex)-1)] = matrix(as.numeric(unlist(gtex[, 3:11690])),ncol = 11688)

# pre-calculate all pairwise correlations
# sample from this
all.RR = cor(t(gtex[,3:11690]))

# find edges in gtex
Ia = match(ensA, gtex$ensembl)
Ib = match(ensB, gtex$ensembl)

# get co-expression correlation for edges
edges.RR = numeric(length(edges))
Igood = which(!is.na(Ia) & !is.na(Ib))
for (ii in 1:length(Igood)) {
  edges.RR[Igood[ii]] = all.RR[Ia[Igood[ii]], Ib[Igood[ii]]]
}
edges.RR[edges.RR==0] = NA




#print("plotting...")
# plot
# GO
# I.novel = (!edges.in.interactome) & both.entrez==1
# I.known = (edges.in.interactome) & both.entrez==1
# N.novel = sum(I.novel)
# N.known = sum(I.known)
# onts = c("BP","CC","MF")
# df = data.frame(ont=character(6), implied.or.not=character(6),
#                 fraction.gt0=numeric(6), mean=numeric(6),
#                 stringsAsFactors = F)
# for (ii in 1:length(onts)) {
#   df$ont[ii] = onts[ii]
#   df$implied.or.not[ii] = "Implied"
#   df$fraction.gt0[ii] = sum(N.terms[,ii]>0 & I.novel, na.rm=T) / N.novel
#   df$mean[ii] = mean(N.terms[I.novel,ii], na.rm=T)
#   
#   df$ont[ii+3] = onts[ii]
#   df$implied.or.not[ii+3] = "Interactome"
#   df$fraction.gt0[ii+3] = sum(N.terms[,ii]>0 & I.known, na.rm=T) / N.known
#   df$mean[ii+3] = mean(N.terms[I.known,ii], na.rm=T)
# }
# ggplot(df, aes(x=ont, y=fraction.gt0)) + 
#   geom_bar(stat = "identity",position = "dodge",aes(fill = implied.or.not))
# ggsave("../figures/impliedEdges_fraction1GO.png")
# 
# ggplot(df, aes(x=ont, y=mean)) + 
#   geom_bar(stat = "identity",position = "dodge",aes(fill = implied.or.not))
# ggsave("../figures/impliedEdges_avgsharedGO.png")
# 
# 
# # Co-expression
# df = data.frame(R=edges.RR, implied.or.not=!edges.in.interactome, stringsAsFactors = F)
# ggplot(df, aes(x=R, fill=implied.or.not)) + geom_density(alpha = 0.3)
# ggsave("../figures/impliedEdges_coexpress.png")



print("write data...")
# write to file
df2write = data.frame(edges=edges, 
                      RR=edges.RR, 
                      both.entrez, 
                      in.interactome=edges.in.interactome,
                      N.terms.BP=N.terms$BP,
                      N.terms.CC=N.terms$CC,
                      N.terms.MF=N.terms$MF, stringsAsFactors = F)
fn = paste("../data/impliedEdges_",
           noise_type,"_",noise_mag,".txt", sep="")
write_tsv(df2write, fn)

