
require(ggplot2)
source("functions.R")

#fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/data/fig2B_v05.Rda"
#load(fn) # sim, Ji
fn = "../data/clusters_full_netw_walktrap.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[!Ji$algorithm=="hierarchical",]
Ji = Ji[Ji$noise_mag %in% c(0,0.01,0.02,0.05,0.1,0.15,0.25,0.5,1),]

Ji$measure = Ji$algorithm
Ji$algorithm[Ji$algorithm=="co_mcl"] = "CO+MCL"
Ji$algorithm[Ji$algorithm=="co"] = "CO"
Ji$algorithm[Ji$algorithm=="mcl"] = "MCL"
Ji$algorithm[Ji$algorithm=="pam"] = "k-Med"
Ji$algorithm[Ji$algorithm=="walk"] = "walktrap"

# make cluster.size, remove all clusters with size<3
Ji$cluster.size = sapply((sapply(Ji$cluster, strsplit, ";")), length)
Ji = Ji[Ji$cluster.size>2,]

# make size.factor
Ji$noise_mag = as.numeric(Ji$noise_mag)
Ji$size.factor = character(nrow(Ji))
Ji$size.factor[Ji$cluster.size<=3] = "N<=3"
Ji$size.factor[Ji$cluster.size<=6 & Ji$cluster.size>3] = "3<N<=6"
Ji$size.factor[Ji$cluster.size<=12 & Ji$cluster.size>6] = "6<N<=12"
Ji$size.factor[Ji$cluster.size>12] = "N>12"
Ji$size.factor = factor(Ji$size.factor, levels= c("N<=3", "3<N<=6", "6<N<=12", "N>12"))

# make size.prctile
Ji$size.prctile = numeric(nrow(Ji))
for (ii in 1:length(unique(Ji$noise_mag))) {
  for (jj in 1:length(unique(Ji$algorithm))) {
    for (kk in 1:length(unique(Ji$iter))) {
      I = Ji$noise_mag==sort(unique(Ji$noise_mag))[ii] & 
        Ji$algorithm==sort(unique(Ji$algorithm))[jj] & 
        Ji$iter==sort(unique(Ji$iter))[kk]
      Ji$size.prctile[I] = rank(Ji$cluster.size[I]) / sum(I)
      #Ji$size.prctile[I] = quantile(Ji$cluster.size[I], seq(from=0, to=1, length=sum(I)))
    }
  }
}

# make noise.factor
#Ji$noise.factor = factor("fpr=0", levels = c("fpr=0", "fpr<=0.02","0.05<fpr<=0.1","0.15<fpr<=0.25","0.50>=fpr"))
Ji$noise.factor[Ji$noise_mag == 0] = ("fpr=0")
Ji$noise.factor[Ji$noise_mag %in% c(0.01, 0.02)] = ("fpr<=0.02")
Ji$noise.factor[Ji$noise_mag %in% c(0.05, 0.1)] = ("0.05<fpr<=0.1")
Ji$noise.factor[Ji$noise_mag %in% c(0.15, 0.25)] = ("0.15<fpr<=0.25")
Ji$noise.factor[Ji$noise_mag %in% c(0.5, 1)] = ("0.50>=fpr")
#levels(Ji$noise.factor) = c("fpr=0", "fpr<=0.02","0.05<fpr<=0.1","0.15<fpr<=0.25","0.50>=fpr")


fn = "../data/clusters_Ai_vs_fdr_df_z.Rda"
load(fn) # df.z


# get averages for figures
nn = 10^3
dm = data.frame(noise_mag = numeric(nn), 
                Ji1 = numeric(nn),
                Ji2 = numeric(nn),
                algorithm = character(nn),
                size.factor = character(nn),
                stringsAsFactors = F)
cc = 0
for (ii in 1:length(unique(Ji$noise_mag))) {
  this.noise = sort(unique(Ji$noise_mag))[ii]
  for (jj in 1:length(unique(Ji$algorithm))) {
    this.algorithm = sort(unique(Ji$algorithm))[jj]
    for (kk in 1:length(unique(Ji$size.factor))) {
      this.size = sort(unique(Ji$size.factor))[kk]
      I = Ji$noise_mag==this.noise & 
        Ji$algorithm==this.algorithm & 
        Ji$size.factor==this.size
      if (sum(I)<10) next
      cc = cc+1
      dm[cc,] = c(this.noise, mean(Ji$Ji1[I], na.rm=T), mean(Ji$Ji2[I], na.rm=T), this.algorithm, this.size)
    }
  }
}
dm = dm[1:cc,]
dm = dm[!is.na(dm$size.factor),]
dm$noise_mag = as.numeric(dm$noise_mag)
dm$Ji1 = as.numeric(dm$Ji1)
dm$Ji2 = as.numeric(dm$Ji2)
tmp = sort(unique(Ji$size.factor))
dm$size.factor = tmp[as.numeric(dm$size.factor)]
dm = dm[!dm$algorithm=="hierarchical",]


# get null for panel B
fn = "../data/clusters_full_netw_walktrap_null.txt"
Ji.null = as.data.frame(read_tsv(fn))
names(Ji.null) = c("cluster.iter", "scramble.iter", "algorithm", "cluster", "Ji1")
Ji.null$cluster.size = sapply((sapply(Ji.null$cluster, strsplit, ";")), length)
Ji.null = Ji.null[Ji.null$cluster.size>2, ]
Ji.null$algorithm[Ji.null$algorithm=="co_mcl"] = "CO+MCL"
Ji.null$algorithm[Ji.null$algorithm=="co"] = "CO"
Ji.null$algorithm[Ji.null$algorithm=="mcl"] = "MCL"
Ji.null$algorithm[Ji.null$algorithm=="pam"] = "k-Med"
Ji.null$algorithm[Ji.null$algorithm=="walk"] = "walktrap"
# make summary
unqsize = sort(unique(Ji.null$cluster.size))
unqalgs = unique(Ji.null$algorithm)
nn = 10^4
null.summ = data.frame(cluster.size = integer(nn),
                       algorithm = character(nn),
                       avg = numeric(nn),
                       avg.95 = numeric(nn),
                       avg.05 = numeric(nn),
                       nn = numeric(nn), stringsAsFactors = F)
cc = 0
for (jj in 1:length(unqalgs)) {
  for (ii in 1:length(unqsize)) {
    cc = cc+1
    I = Ji.null$cluster.size==unqsize[ii] & Ji.null$algorithm==unqalgs[jj]
    null.summ$cluster.size[cc] = unqsize[ii]
    null.summ$algorithm[cc] =unqalgs[jj]
    null.summ$avg[cc] = mean(Ji.null$Ji1[I], na.rm = T)
    null.summ$avg.95[cc] = quantile(Ji.null$Ji1[I], .975)
    null.summ$avg.05[cc] = quantile(Ji.null$Ji1[I], .025)
    null.summ$nn[cc] = sum(I)
  }
  I = null.summ$algorithm == unqalgs[jj]
  null.summ$avg[I] = predict(loess(avg ~ cluster.size, data = null.summ[I,]), newdata = null.summ[I,])
  null.summ$avg.95[I] = predict(loess(avg.95 ~ cluster.size, data = null.summ[I,]), newdata = null.summ[I,])
  null.summ$avg.05[I] = predict(loess(avg.05 ~ cluster.size, data = null.summ[I,]), newdata = null.summ[I,])
}
null.summ = null.summ[1:cc, ]

# 4B. Cluster size vs Ji
Ji$noise.factor = as.factor(Ji$noise.factor)
Ji$noise.factor = ordered(Ji$noise.factor, 
                          levels = c("0.50>=fpr", "0.15<fpr<=0.25","0.05<fpr<=0.1","fpr<=0.02", "fpr=0"))
Ji$alpha = I(0.005 + 0.09 * Ji$cluster.size / 100)
Ji$alpha[Ji$alpha>0.5] = 0.5
Ji$alpha[Ji$Ji2==1] = Ji$alpha[Ji$Ji2==1] / 3
ggplot(data = Ji[Ji$iter>1,], aes(x=log10(cluster.size), y=Ji2, color=noise.factor, alpha=alpha)) + 
  geom_line(data = null.summ, inherit.aes = F, aes(x=log10(cluster.size), y=avg), color="black") + 
  geom_line(data = null.summ, inherit.aes = F, aes(x=log10(cluster.size), y=avg.95), color="black", linetype="dashed") + 
  geom_line(data = null.summ, inherit.aes = F, aes(x=log10(cluster.size), y=avg.05), color="black", linetype="dashed") + 
  geom_point() + 
  facet_grid(~algorithm) +
  geom_smooth(method="lm") + 
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
  scale_x_continuous(breaks=log10(c(3,10,25,50,100,200)), labels = c(3,10,25,50,100,200)) +
  coord_cartesian(ylim = c(-.001,1.001), xlim=c(log10(2.6),log10(200))) + theme_bw() + 
  scale_color_viridis(discrete = T) + 
  theme(legend.position = "none") +
  ylab("Similarity to iter=1 (Ji)") + 
  xlab("Cluster size") 
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4C_v04.png"
ggsave(fn,width=10, height=3)

ggplot() + 
  geom_point(data = Ji.null, aes(x=log10(cluster.size), y=Ji1), color="red", alpha = .01) + 
  facet_grid(~algorithm)

# 4D. Consensus adjacency matrix examples
load("/Users/gregstacey/Academics/Foster/ClusterExplore/data/clusters_Ai_vs_fdr.Rda")
I0 = clusters$algorithm=="mcl" & clusters$noise_mag==0.15
clusters = clusters[I0,]
I1 = which(clusters$iter==1)

clust.good = clusters$cluster[166]
clust.bad =  clusters$cluster[78]

df.adj.good = consensus.adjmat(clust.good, clusters)
df.adj.bad = consensus.adjmat(clust.bad, clusters)
df.adj.good$iter.name = factor(paste("iter=", df.adj.good$iter, sep=""))
df.adj.bad$iter.name = factor(paste("iter=", df.adj.bad$iter, sep=""))
df.adj.good$iter.name = factor(df.adj.good$iter.name, levels=paste("iter=", 1:11, sep=""))
df.adj.bad$iter.name = factor(df.adj.bad$iter.name, levels=paste("iter=", 1:11, sep=""))

df.adj.good = df.adj.good[df.adj.good$iter %in% c(1,5,10),]
ggplot(df.adj.good, aes(prots, variable)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  facet_grid(~iter.name) + theme(legend.position="none") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4D_v03.png"
ggsave(file=fn, width=10*3/5 * .8, height=2.3 * .8)

df.adj.bad = df.adj.bad[df.adj.bad$iter %in% c(1,5,10),]
ggplot(df.adj.bad, aes(prots, variable)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  facet_grid(~iter.name) + theme(legend.position="none") +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank(),axis.text.y=element_blank(),legend.position = "none",
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_4E_v03.png"
ggsave(file=fn, width=10*3/5 * .8, height=2.3 * .8)



# enrichment of good and bad clusters
# read GO
fn = "../data/gene2go"
gene2go = read_tsv(fn)
gene2go = gene2go[!gene2go$Evidence %in% "IEA", ] # filter by evidence

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

# load clusters
fn = "../data/clusters_full_netw.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[Ji$iter==1,]

# protein universe
this.background = unique(unlist(strsplit(Ji$cluster, ";")))
this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
n.universe = length(this.background.entrez)

# get GO associated with this universe
this.gene2go = gene2go[gene2go$GeneID %in% this.background.entrez,]

# loop through ontologies
onts = unique(this.gene2go$Category)

clusts = list()
clusts[[1]] = clust.good
clusts[[2]] = clust.bad
df = data.frame(clust = numeric(10^3),
                ont = character(10^3),
                p.or.q = character(10^3), 
                go.terms = character(10^3), stringsAsFactors = F)
cc = 0
for (ii in 1:length(clusts)) {
  for (oo in 1:length(onts)) {
    ont = onts[oo]
    
    # get GO associated with this universe and ontology
    this.gene2go.ont = this.gene2go[this.gene2go$Category == ont,]
    
    # remove redundant entries from this.gene2go, and filter to >20 & <200 proteins
    this.gene2go.ont = distinct(this.gene2go.ont[,c("GeneID", "GO_ID")])
    x = table(this.gene2go.ont[,c("GeneID", "GO_ID")])
    n.times = colSums(x)
    filtered.terms = colnames(x)[n.times>10 & n.times<200]
    this.gene2go.ont = this.gene2go.ont[this.gene2go.ont$GO_ID%in%filtered.terms,]
    this.goterms = unique(this.gene2go.ont$GO_ID)
    
    # filter to terms annotated 
    nn = length(this.goterms)
    
    # set
    this.ids = unlist(strsplit(clusts[[ii]], ";"))
    this.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
    n.set = length(this.entrez)
    if (length(this.entrez)<2) {
      print("    cluster size <= 2...")
      next
    }
    
    sumTable = data.frame(GOBPID = character(nn), 
                          Pvalue = numeric(nn),
                          n.u = numeric(nn),n.u.hits = numeric(nn),n.s = numeric(nn),n.s.hits = numeric(nn),
                          stringsAsFactors = F)
    for (mm in 1:nn) {
      go = this.goterms[mm]
      sumTable$GOBPID[mm] = go
      
      # universe hits
      n.u.hits = sum(this.gene2go.ont$GeneID%in%this.background.entrez & this.gene2go.ont$GO_ID %in% go)
      
      # set hits
      n.s.hits = sum(this.gene2go.ont$GeneID %in% this.entrez & this.gene2go.ont$GO_ID %in% go)
      
      sumTable$Pvalue[mm] = phyper(n.s.hits-1, n.u.hits, 
                                   n.universe-n.u.hits, n.set, lower.tail=FALSE)
      sumTable$n.u[mm] = n.universe
      sumTable$n.u.hits[mm] = n.u.hits
      sumTable$n.s[mm] = n.set
      sumTable$n.s.hits[mm] = n.s.hits
    }
    sumTable[is.na(sumTable)] = 1
    sumTable$qvalue = p.adjust(sumTable$Pvalue)
    
    # p-significant
    cc = cc+1
    I.p = sumTable$Pvalue <= 0.05
    df$clust[cc] = ii
    df$ont[cc] = ont
    df$p.or.q[cc] = "p"
    df$go.terms[cc] = paste(sumTable$GOBPID[I.p], collapse=";")
    
    # q-significant
    cc = cc+1
    I.q = sumTable$qvalue <= 0.05
    df$clust[cc] = ii
    df$ont[cc] = ont
    df$p.or.q[cc] = "q"
    df$go.terms[cc] = paste(sumTable$GOBPID[I.q], collapse=";")
    
  }
}
df = df[1:cc,]

# how many of the good cluster proteins are ribosomal?
#fn = "../data/uniprot-human-ribosome.txt"
#hurib = as.data.frame(read_tsv(fn))
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))


this.ids = unlist(strsplit(clust.good, ";"))
I = grep("ibosom", corum$ComplexName)
for (ii in 1:length(I)) {
  prots = unlist(strsplit(corum$`subunits(UniProt IDs)`[I[ii]], ";"))
  nn = sum(this.ids %in% prots)
  name = corum$ComplexName[I[ii]]
  print(paste(nn, "/", length(this.ids), "in", name))
}

this.ids = unlist(strsplit(clust.bad, ";"))
JJ = numeric(nrow(corum))
nn = numeric(nrow(corum))
for (ii in 1:nrow(corum)) {
  prots = unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";"))
  nn[ii] = sum(this.ids %in% prots)
  JJ[ii] = sum(this.ids %in% prots) / length(unique(c(nn, prots)))
  print(paste(ii, nn[ii]))
}



allComplexes$Protein.IDs <- character(length(allComplexes$subunits.UniProt.IDs)) # preallocate column
for (ii in 1:nrow(allComplexes)) { # go through every protein you want to match
  these.prots = unlist(strsplit(as.character(allComplexes$subunits.UniProt.IDs[ii]),";")) # split into single IDs
  human.matches = character(length(these.prots)) # preallocate mouse matches
  for (jj in 1:length(these.prots)) { # go through every single ID of that protein
    I = which(H.sapiens.M.musculus$Human %in% these.prots[jj])
    if (length(I)==1) { # if you found exactly ONE human-to-mouse ortholog...
      human.matches[jj] = H.sapiens.M.musculus$Protein.IDs[I] # ... store the mouse ID in human matches
    }
  }
  # concatenate mouse matches
  human.matches = human.matches[ !human.matches %in% ""] # remove NAs
  human.matches = paste(human.matches, collapse=";") # collapse into semicolon separated string
  
  # store in condition file
  allComplexes$Protein.IDs[ii] = human.matches
}

