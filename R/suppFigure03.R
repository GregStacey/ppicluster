
source("functions.R")
require(dils)

# matlab results (chromatograms + network, three algorithms)
fn = "../data/clusters_full_netw_walktrap.txt"
Ji = as.data.frame(read_tsv(fn))
Ji = Ji[!Ji$algorithm=="hierarchical",]
Ji = Ji[Ji$noise_mag %in% c(0,0.01,0.02,0.05,0.1,0.15,0.25,0.5,1),]


# Count the rearrangements. "Simple" analysis
sf = "../data/dfchange_02.Rda"
if (F) {
  load(sf)
} else {
  df.change = data.frame(noise = numeric(10^4),
                         n.shuffled.interactions = numeric(10^4),
                         algorithm = character(10^4),
                         n.new.clusters = numeric(10^4), 
                         n.lost.clusters = numeric(10^4),
                         n.clusters = numeric(10^4),
                         avg.gained.nodes = numeric(10^4),
                         avg.lost.nodes = numeric(10^4), 
                         total.edges0 = numeric(10^4),
                         n.rearranged.edges = numeric(10^4),
                         stringsAsFactors = F)
  
  unqiter = unique(Ji$iter)
  unqnoise = sort(unique(Ji$noise_mag))
  unqalgs = unique(Ji$algorithm)
  cc = 0
  
  for (kk in 1:length(unqiter)) {
    print(paste("iter =", unqiter[kk]))
    for (mm in 1:length(unqalgs)) {
      set0 = Ji$cluster[Ji$noise_mag==0 & Ji$algorithm%in%unqalgs[mm] & Ji$iter==unqiter[kk]]
      print(unqalgs[mm])
      for (uu in 1:length(unqnoise)) {
        print(unqnoise[uu])
        set1 = Ji$cluster[Ji$noise_mag==unqnoise[uu] & Ji$algorithm%in%unqalgs[mm] & Ji$iter==unqiter[kk]]
        
        # scan set1 for n.gained
        n.gained.nodes = numeric(length(set1))
        n.set1 = numeric(length(set1))
        n.edges1 = numeric(length(set1))
        for (ii in 1:length(set1)) {
          # find best-matching cluster
          cluster1 = unlist(strsplit(set1[ii], ";"))
          JJ = numeric(length(set0))
          for (jj in 1:length(set0)) {
            cluster0 = unlist(strsplit(set0[jj], ";"))
            nn.jacc = length(unique(c(cluster0, cluster1)))
            JJ[jj] = length(intersect(cluster0, cluster1)) / nn.jacc 
          }
          I.max = which.max(JJ)[1]
          cluster0 = unlist(strsplit(set0[I.max], ";"))
          n.gained.nodes[ii] = sum(!cluster1 %in% cluster0)
          n.set1[ii] = length(cluster1)
          n.edges1[ii] = length(cluster1) * (length(cluster1)-1) / 2
        }
        
        # scan set0 for n.lost.clusters
        n.lost.nodes = numeric(length(set0))
        n.set0 = numeric(length(set0))
        n.edges0 = numeric(length(set0))
        n.shuffle.edges0 = numeric(length(set0))
        for (ii in 1:length(set0)) {
          # find best-matching cluster
          cluster0 = unlist(strsplit(set0[ii], ";"))
          JJ = numeric(length(set1))
          for (jj in 1:length(set1)) {
            cluster1 = unlist(strsplit(set1[jj], ";"))
            nn.jacc = length(unique(c(cluster1, cluster0)))
            JJ[jj] = length(intersect(cluster1, cluster0)) / nn.jacc 
          }
          I.max = which.max(JJ)[1]
          cluster1 = unlist(strsplit(set1[I.max], ";"))
          n.lost.nodes[ii] = sum(!cluster0 %in% cluster1)
          n.set0[ii] = length(cluster0)
          n.edges0[ii] = length(cluster0) * (length(cluster0)-1) / 2
        }
        
        cc = cc+1
        df.change$noise[cc] = unqnoise[uu]
        df.change$n.shuffled.interactions[cc] = round(39563 * unqnoise[uu])
        df.change$algorithm[cc] = unqalgs[mm]
        df.change$n.lost.clusters[cc] = sum((n.lost.nodes-n.set0)==0)
        df.change$n.new.clusters[cc] = sum((n.gained.nodes-n.set1)==0)
        df.change$n.clusters[cc] = length(set1)
        df.change$avg.gained.nodes[cc] = mean(n.gained.nodes)
        df.change$avg.lost.nodes[cc] = mean(n.lost.nodes)
        df.change$total.edges0[cc] = sum(n.edges0)
        
        # calculate rearranged edges as a difference between adjacency matrices
        # edgelist, cluster0
        nn = sum(n.edges0)
        elist0 = data.frame(protA = character(nn), protB = character(nn), weights = numeric(nn)+1, stringsAsFactors = F)
        cc2 = 0
        for (ii in 1:length(set0)) {
          cluster0 = unlist(strsplit(set0[ii], ";"))
          if (length(cluster0)<3) next
          ab = combn(cluster0, 2)
          I = (cc2+1) : (cc2+ncol(ab))
          elist0$protA[I] = ab[1,]
          elist0$protB[I] = ab[2,]
          cc2 = cc2+ncol(ab)
        }
        # edgelist, cluster1
        nn = sum(n.edges1)
        elist1 = data.frame(protA = character(nn), protB = character(nn), weights = numeric(nn)+1, stringsAsFactors = F)
        cc2 = 0
        for (ii in 1:length(set1)) {
          cluster1 = unlist(strsplit(set1[ii], ";"))
          if (length(cluster1)<3) next
          ab = combn(cluster1, 2)
          I = (cc2+1) : (cc2+ncol(ab))
          elist1$protA[I] = ab[1,]
          elist1$protB[I] = ab[2,]
          cc2 = cc2+ncol(ab)
        }
        # make common indices, sort
        allprots = unique(c(elist0$protA, elist0$protB, elist1$protA, elist1$protB))
        ia0 = match(elist0$protA, allprots)
        ib0 = match(elist0$protB, allprots)
        elist0.num = cbind(ia0, ib0)
        elist0.num = t(apply(elist0.num, 1, sort))
        ia1 = match(elist1$protA, allprots)
        ib1 = match(elist1$protB, allprots)
        elist1.num = cbind(ia1, ib1)
        elist1.num = t(apply(elist1.num, 1, sort))
        if (nrow(elist1.num)<2) elist1.num = cbind(0,0)
        # count differences
        e0 = elist0.num[,1] + elist0.num[,2]*10^4
        e1 = elist1.num[,1] + elist1.num[,2]*10^4
        # add occurrence number to count duplicates
        n0 = as.data.frame(e0) %>% group_by(e0) %>% mutate(Index=1:n())
        n1 = as.data.frame(e1) %>% group_by(e1) %>% mutate(Index=1:n())
        e0 = e0 + n0$Index*10^8
        e1 = e1 + n1$Index*10^8
        
        df.change$n.rearranged.edges[cc] = length(unique(c(e0, e1))) - length(intersect(e0, e1))
      }
    }
  }
  df.change = df.change[1:cc,]
  
  save(df.change, file = sf)
}



# number of clusters
ggplot(df.change, aes(x=noise, y=n.clusters, color=algorithm)) + geom_line() +
  scale_colour_grey() + xlab("Interactome FPR") + ylab("Number of clusters") +  
  theme_bw()
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup3E1_v02.png"
ggsave(fn,width=4.5, height=3)

# avg gained and lost proteins
df = melt(df.change, id.vars=c("noise", "algorithm"), measure.vars = c("avg.gained.nodes", "avg.lost.nodes"))
df$variable = as.character(df$variable)
df[df=="avg.gained.nodes"] = "Gained"
df[df=="avg.lost.nodes"] = "Lost"
ggplot(df, aes(x=noise, y=value, color=algorithm)) + geom_line() + 
  facet_wrap(~variable, scales="free_y") + theme_bw() +
  scale_colour_grey() + xlab("Interactome FPR") + 
  ylab("Number of gained or lost proteins\nper cluster (avg.)") +
  theme(legend.position = "none")
fn = "/Users/gregstacey/Academics/Foster/Manuscripts/ClusterExplore/figures/fig_sup3E2_v02.png"
ggsave(fn,width=6.5, height=3)

# number of gained and lost clusters
df = melt(df.change, id.vars=c("noise", "algorithm"), measure.vars = c("n.new.clusters", "n.lost.clusters"))
df$variable = as.character(df$variable)
df[df=="n.new.clusters"] = "Gained"
df[df=="n.lost.clusters"] = "Lost"
ggplot(df, aes(x=noise, y=value, color=algorithm)) + geom_line() + 
  facet_wrap(~variable, scales="free_y") +
  scale_colour_grey() + xlab("Interactome FPR") + ylab("Number of gained or lost clusters")




# number of shuffled interactions vs number of rearranged edges
# N_rearranged = N.lost + N.gained
