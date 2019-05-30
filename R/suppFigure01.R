
# A - network connection matrices

source("functions.R")

# read corum
fn = "../data/allComplexes.txt"
corum = as.data.frame(read_tsv(fn))
corum = corum[corum$Organism=="Human",]

# turn corum into a list of protein pairs
ints = list()

ints[[1]] = data.frame(protA = character(10^6),
                        protB = character(10^6), stringsAsFactors = F)
cc = 0
for (ii in 1:nrow(corum)) {
  print(ii)
  prots = sort(unlist(strsplit(corum$`subunits(UniProt IDs)`[ii], ";")))
  if (length(prots)<2) next
  pairs.prots = t(combn(prots, 2))
  
  I = (cc+1) : (cc+nrow(pairs.prots))
  ints[[1]]$protA[I] = pairs.prots[,1]
  ints[[1]]$protB[I] = pairs.prots[,2]
  cc = cc+length(I)
}
ints[[1]] = ints[[1]][1:cc,]
ints[[1]] = distinct(ints[[1]])


# read chem-protein and facebook networks
fns = c("", "../data/ChCh-Miner_durgbank-chem-chem.tsv", "../data/email-Eu-core.txt")
seps = c("", "\t", " ")
for (ii in 2:length(fns)) {
  ints[[ii]] = as.data.frame(read_delim(fns[ii], delim=seps[ii], quote = ""))
  colnames(ints[[ii]]) = c("protA","protB")
}
ints[[3]] = ints[[3]] + 1


# convert to adjacency matrices
df = list()
for (ii in 1:3) {
  print(ii)
  prots = sort(unique(c(ints[[ii]]$protA, ints[[ii]]$protB)))
  conmat = matrix(nrow=length(prots), ncol=length(prots))
  ia = match(ints[[ii]]$protA, prots)
  ib = match(ints[[ii]]$protB, prots)
  for (jj in 1:length(ia)) {
    conmat[ia[jj], ib[jj]] = 1
    conmat[ib[jj], ia[jj]] = 1
  }
  
  conmat.df = as.data.frame(conmat)
  names(conmat.df) = as.character(1:length(prots))
  conmat.df$proteins = as.character(1:length(prots))
  df[[ii]] = meltmat(conmat.df, id.vars = c("proteins"))
  df[[ii]]$value[is.na(df[[ii]]$value)] = 0
}


ii = 1
ggplot(df[[ii]], aes(x, y)) +
  geom_tile(aes(fill = value), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  theme_bw() + blank_theme + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

