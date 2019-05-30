source("functions.R")


if (F){
  fn = "../data/clusters_full_netw_wEnr2.txt"
  Ji = as.data.frame(read_tsv(fn))
  Ji = Ji[Ji$iter==1,]
  
  fn = "../data/clusters_full_netw2.txt"
  Ji0 = as.data.frame(read_tsv(fn))
  Ji0 = Ji0[Ji0$iter==1,]
  
  Ji = cbind(Ji0[,c("iter","noise_mag","algorithm","cluster","Ji1","Ji2")],
             Ji[,c("bp.enriched.p","cc.enriched.p","mf.enriched.p",
                   "bp.enriched.q","cc.enriched.q","mf.enriched.q")])
} else {
  
  # load clusters
  fn = "../data/clusters_full_netw.txt"
  Ji = as.data.frame(read_tsv(fn))
  Ji = Ji[Ji$iter==1,]
  
  # save unadulterated Ji for later
  Ji0 = Ji
  
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
  
  # calculate GO enrichment
  Ji$bp.enriched.p = character(nrow(Ji))
  Ji$cc.enriched.p = character(nrow(Ji))
  Ji$mf.enriched.p = character(nrow(Ji))
  Ji$bp.enriched.q = character(nrow(Ji))
  Ji$cc.enriched.q = character(nrow(Ji))
  Ji$mf.enriched.q = character(nrow(Ji))
  ont.columns.p = c("bp.enriched.p","cc.enriched.p","mf.enriched.p")
  I.ont.columns.p = match(ont.columns.p,names(Ji))
  ont.columns.q = c("bp.enriched.q","cc.enriched.q","mf.enriched.q")
  I.ont.columns.q = match(ont.columns.q,names(Ji))
  unqnoise = unique(Ji$noise_mag)
  unqalgs = unique(Ji$algorithm)
  
  for (ii in 1:nrow(Ji)) {
    I.this.set = which(Ji$noise_mag==Ji$noise_mag[ii] & Ji$algorithm==Ji$algorithm[ii])
    
    # protein universe
    this.background = unique(unlist(strsplit(Ji$cluster[I.this.set], ";")))
    this.background.entrez = unique(ens_ent$entrez[ens_ent$uniprot %in% this.background])
    n.universe = length(this.background.entrez)
    
    # get GO associated with this universe
    this.gene2go = gene2go[gene2go$GeneID %in% this.background.entrez,]
    
    # loop through ontologies
    onts = unique(this.gene2go$Category)
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
      this.ids = unlist(strsplit(Ji$cluster[ii], ";"))
      this.entrez = ens_ent$entrez[ens_ent$uniprot %in% this.ids]
      n.set = length(this.entrez)
      if (length(this.entrez)<2) {
        print("    cluster size <= 2...")
        next
      }
      
      print(paste(ii, Ji$algorithm[ii], Ji$noise_mag[ii]))
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
      I.p = sumTable$Pvalue <= 0.05
      Ji[ii, I.ont.columns.p[oo]]= paste(sumTable$GOBPID[I.p], collapse=";")
      
      # q-significant
      I.q = sumTable$qvalue <= 0.05
      Ji[ii, I.ont.columns.q[oo]]= paste(sumTable$GOBPID[I.q], collapse=";")
      
      if (sum(I.p)>0 | sum(I.q)>0) {
        print(paste(" !!  N.sig.p = ", sum(I.p), ",  N.sig.q = ", sum(I.q)))
      }
      
    }
  }
  
  # save Ji before binding
  fn = "../data/clusters_full_netw_wEnr.Rda"
  save(Ji, file=fn)
  
  Ji = cbind(Ji0[,c("iter","noise_mag","algorithm","cluster","Ji1","Ji2")],
             Ji[,c("bp.enriched.p","cc.enriched.p","mf.enriched.p",
                   "bp.enriched.q","cc.enriched.q","mf.enriched.q")])
  
  # write finally
  fn = "../data/clusters_full_netw_wEnr.txt"
  write_tsv(Ji, path=fn)
}



Ji$nsig.p = numeric(nrow(Ji))
Ji$nsig.q = numeric(nrow(Ji))
for (ii in 1:nrow(Ji)) {
  x1 = length(unlist(strsplit(Ji$bp.enriched.p[ii], ";")))
  x2 = length(unlist(strsplit(Ji$mf.enriched.p[ii], ";")))
  x3 = length(unlist(strsplit(Ji$cc.enriched.p[ii], ";")))
  x4 = length(unlist(strsplit(Ji$bp.enriched.q[ii], ";")))
  x5 = length(unlist(strsplit(Ji$mf.enriched.q[ii], ";")))
  x6 = length(unlist(strsplit(Ji$cc.enriched.q[ii], ";")))
  Ji$nsig.p[ii] = x1+x2+x3
  Ji$nsig.q[ii] = x4+x5+x6
}

ggplot(Ji, aes(x=noise_mag, y=nsig.p)) + geom_jitter(alpha=.05)