
require(grex)
require(ggplot2)

fnsave = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters_coexp.Rda"
load(fnsave)


# load clusters
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/clusters.txt"
clusters = read.csv2(fn, sep="\t", quote="", stringsAsFactors = F)
clusters = clusters[,!names(clusters) %in% "X"]
clusters$noise_mag = as.numeric(clusters$noise_mag)

# add real CORUM complexes as positive control
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/allComplexes.txt"
corum = read.table(fn, sep="\t", stringsAsFactors = F, quote="", header = T)
corum2add = data.frame(data_type=character(nrow(corum)),
                       noise_type=character(nrow(corum)),
                       noise_mag=numeric(nrow(corum)),
                       algorithm=character(nrow(corum)),
                       cluster=corum$subunits.UniProt.IDs.,
                       stringsAsFactors = F)
corum2add$data_type = "corum_real"
corum2add$noise_type = "-1"
corum2add$noise_mag = "-1"
corum2add$algorithm = "-1"
clusters = rbind(clusters, corum2add)


# unique IDs
unqprots = unique(unlist(strsplit(clusters$cluster, ";")))



# Co-expression

# read gtex data
fn = "/Users/Mercy/Academics/Foster/ClusterExplore/data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz"
gtex = read.table(fn, sep="\t", quote="", stringsAsFactors = F, skip = 3)
gtex$ensembl = sapply(strsplit(gtex$V1, split=".", fixed=T), "[", 1)

# map uniprot_clusters
data(gtexv7)
df.idmap = grex(gtexv7)
ens_mappable  = df.idmap$ensembl_id[df.idmap$uniprot_id %in% unqprots]
gtex = gtex[gtex$ensembl %in% ens_mappable,]
clusters$ensembl = character(nrow(clusters))
clusters$size.ensembl = numeric(nrow(clusters))
for (ii in 1:nrow(clusters)){
  this.cluster = unlist(strsplit(clusters$cluster[ii], ";", fixed=T))
  I = df.idmap$uniprot_id %in% this.cluster
  clusters$size.ensembl[ii] = sum(I)
  if (sum(I)<2) next
  clusters$ensembl[ii] = paste(df.idmap$ensembl_id[I], collapse=";")
}
# dumb correction: make gtex columns numeric
gtex[,3:(ncol(gtex)-1)] = matrix(as.numeric(unlist(gtex[, 3:11690])),ncol = 11688)

# pre-calculate all pairwise correlations
# sample from this
all.RR = cor(t(gtex[,3:11690]))

# calculate pairwise correlation
I = which(!clusters$ensembl %in% "")
clusters$R.avg.pairwise = numeric(nrow(clusters))
clusters$R.zscore = numeric(nrow(clusters))
clusters$R.rand = numeric(nrow(clusters))
clusters$R.avg.pairwise[clusters$ensembl %in% ""] = NA
clusters$R.zscore[clusters$ensembl %in% ""] = NA
clusters$R.rand[clusters$ensembl %in% ""] = NA
for (ii in 1:length(I)){
  print(ii)
  this.cluster = unlist(strsplit(clusters$ensembl[I[ii]], ";", fixed=T))
  Igtex = gtex$ensembl %in% this.cluster
  RR = all.RR[Igtex, Igtex]
  RR[diag(nrow(RR))==1] = NA
  this.R = mean(as.vector(RR), na.rm=T)
  
  # significance
  iterMax = 100
  rand.RR = numeric(iterMax)
  for (iter in 1:iterMax) {
    fake.I = sample(nrow(gtex), sum(Igtex))
    fake.RR = all.RR[fake.I, fake.I]
    fake.RR[diag(nrow(fake.RR))==1] = NA
    rand.RR[iter] = mean(as.vector(fake.RR), na.rm=T)
  }
  clusters$R.rand[I[ii]] = mean(rand.RR, na.rm=T)
  
  clusters$R.avg.pairwise[I[ii]] = this.R
  clusters$R.zscore[I[ii]] = (this.R - mean(rand.RR, na.rm=T)) / sd(rand.RR, na.rm=T)
}
clusters$noise_mag.factor = as.factor(clusters$noise_mag)




I = !is.na(clusters$R.avg.pairwise) & clusters$data_type %in% "corum"
ggplot(clusters[I,], aes(fill=noise_mag.factor, y=R.avg.pairwise)) + geom_boxplot() + 
  facet_grid(algorithm ~ noise_type)
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/coexp_corum_box_R.png")
ggplot(clusters[I,], aes(fill=noise_mag.factor, y=R.zscore)) + geom_boxplot() + 
  facet_grid(algorithm ~ noise_type)
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/coexp_corum_box_Z.png")
ggplot(clusters[I,], aes(fill=noise_mag.factor, y=log10(size.ensembl))) + geom_boxplot() + 
  facet_grid(algorithm ~ noise_type)
ggsave("/Users/Mercy/Academics/Foster/ClusterExplore/figures/coexp_corum_box_size.png")


I = !is.na(clusters$R.avg.pairwise) & clusters$data_type %in% "corum" & 
  clusters$noise_type %in% "network_remove" & clusters$algorithm %in% "co_mcl"
ggplot(clusters[I,], aes(x=jitter(as.numeric(noise_mag)), y=log10(size.ensembl))) + geom_point(alpha=0.1) + geom_smooth()


I = !is.na(clusters$R.avg.pairwise) & (clusters$data_type %in% "corum" | clusters$data_type %in% "corum_real")
tmp = clusters[I,c("data_type","R.avg.pairwise")]
tmp2 = data.frame(data_type=character(10000), R.avg.pairwise=numeric(10000))
tmp2$data_type = "rand"
tmp2$R.avg.pairwise = clusters$R.rand[sample(nrow(clusters), 10000)]
tmp2 = rbind(tmp, tmp2)
ggplot(tmp2, aes(fill=data_type,x=R.avg.pairwise)) + geom_density(alpha=0.3, bw=0.1)


save(all.RR, clusters, file=fnsave)
