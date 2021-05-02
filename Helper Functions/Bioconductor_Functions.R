
##### ID conversions #####
#Create function to convert between flybase and gene symbols
lookup.id<- function(x){
  suppressMessages(
  AnnotationDbi::select(org.Dm.eg.db, keys=x,columns='SYMBOL',keytype ="FLYBASE")$SYMBOL %>%
    as.character())}



# Function to trim GO terms 
go.trim <- function(go.res, p.cutoff = 0.05) {
  require(GSEABase)
  
  myCollection <- go.res %>%
    filter(over_represented_pvalue <= p.cutoff) %>%
    pull(category) %>%
    GOCollection()
  
  # download go slid database 
  slim <- getOBOCollection('http://current.geneontology.org/ontology/subsets/goslim_drosophila.obo')
  # trim terms 
  goSlim(myCollection, slim, "BP") %>%
    as_tibble(rownames = 'category') %>%
    filter(Count > 0 ) %>%
    return()
}




####### GENE ONTOLOGY #########
library(GO.db)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
# function to convert ids pipeable 
convert.2.fb<- function(x,keytype = 'ENTREZID') {
  AnnotationDbi::mapIds(org.Dm.eg.db,keys = x,column = 'FLYBASE', keytype = keytype, multiVals = 'first') %>%
    as.character()
}

# set up go db conversions
# get table with go id and terms 
# requires txdb object to get gene lengths 
goterm.conversions <- GOTERM %>%
  toTable() 

# get table with genes annotated to go term id 
conversions <- org.Dm.egGO %>% #org.Dm.egGO2ALLEGS
  toTable() %>% 
  mutate(Gene = convert.2.fb(gene_id)) %>% # add column with flybase 
  left_join(goterm.conversions[,2:3]) # merge term tables 



g2cat <- org.Dm.egGO %>%
  toTable() %>%
  mutate(gene_id = convert.2.fb(gene_id)) %>%
  filter(Ontology == 'BP') %>%
  dplyr::select(gene_id, go_id)


get.de <- function(results = results ,p.cutoff = p.cutoff, lfc.cutoff = lfc.cutoff, direction = direction){
  
  if(direction == 'up') {
    de.genes<- results %>%
      filter(padj <= p.cutoff ) %>%
      filter(log2FoldChange > lfc.cutoff) %>%
      pull(Gene) 
    
  } else if (direction == 'down') {
    de.genes<- results %>%
      filter(padj <= p.cutoff ) %>%
      filter(-1 *log2FoldChange > lfc.cutoff ) %>%
      pull(Gene) 
    
  } else if (direction == 'both') {
    de.genes<- results %>%
      filter(padj <= p.cutoff ) %>%
      filter(abs(log2FoldChange) > lfc.cutoff) %>%
      pull(Gene) 
  }
  return(de.genes)
}


gene.ontology2 <- function(results, lfc.cutoff, p.cutoff, correction, direction, fdr.cut) {
  
  de.genes <- get.de(results = results,p.cutoff =p.cutoff , lfc.cutoff = lfc.cutoff,  direction = direction)
  
  # get gene universe
  assayed.genes <- results %>%
    pull(Gene) 
  
  
  #set up list for goseq
  gene.vector=as.integer(assayed.genes %in% de.genes)
  names(gene.vector)=assayed.genes
  
  # Run GO analysis  
  require(goseq)
  pwf = nullp(gene.vector,'dm6','ensGene', plot.fit = F,)
  go.res <-  goseq(pwf,gene2cat = g2cat) %>%
    as_tibble() %>%
    filter(numInCat > 1, numDEInCat > 0)
  
  go.res$padj= p.adjust(go.res$over_represented_pvalue, method= correction) 
  
  go.res <- go.res %>%
    filter(padj <= fdr.cut)
  
  out.genes <- c()
  for (go.term in pull(go.res, category)){
    
    allegs<-get(go.term, org.Dm.egGO2ALLEGS)
    genes = unlist(mget(allegs,org.Dm.egFLYBASE),use.names = FALSE)
    term.cg <- de.genes[which(de.genes %in% genes)] %>%
      lookup.id()
    out.genes <- c(out.genes, toString(term.cg))
  }
  go.res$Symbols <- out.genes
  return(go.res)
}
#############
gene.ontology3 <- function(genes, universe, p.cutoff, correction, fdr.cut) {
  
  de.genes <- genes
  
  # get gene universe
  assayed.genes <-universe

  #set up list for goseq
  gene.vector=as.integer(assayed.genes %in% de.genes)
  names(gene.vector)=assayed.genes
  
  # Run GO analysis  
  require(goseq)
  pwf = nullp(gene.vector,'dm6','ensGene', plot.fit = F,)
  go.res <-  goseq(pwf,gene2cat = g2cat) %>%
    as_tibble() %>%
    filter(numInCat > 1, numDEInCat > 0)
  
  go.res$padj= p.adjust(go.res$over_represented_pvalue, method= correction) 
  
  go.res <- go.res %>%
    filter(padj <= fdr.cut)
  
  out.genes <- c()
  for (go.term in pull(go.res, category)){
    
    allegs<-get(go.term, org.Dm.egGO2ALLEGS)
    genes = unlist(mget(allegs,org.Dm.egFLYBASE),use.names = FALSE)
    term.cg <- de.genes[which(de.genes %in% genes)] %>%
      lookup.id()
    out.genes <- c(out.genes, toString(term.cg))
  }
  go.res$Symbols <- out.genes
  return(go.res)
}
#############

gene.ontology <- function(results, lfc.cutoff = .7, p.cutoff = 0.05, correction = 'BH', separate = TRUE, genes = NULL, universe = NULL, fdr.cut = 0.05){
  if (! is.null(genes)) {
    if (is.null(universe)) stop('Gene universe is required when supplying gene list',call. = F)
    
    go.results = gene.ontology3(genes = genes, universe = universe, p.cutoff = p.cutoff, correction = correction, fdr.cut = fdr.cut)
    
  }
  if (is.null(genes)){
    
  
  if (! separate ) {
    
    go.results = gene.ontology2(results = results, lfc.cutoff = lfc.cutoff, p.cutoff = p.cutoff, correction = correction, direction = 'both', fdr.cut = fdr.cut)
  } 
  else if (separate) {
    go.results.up  = gene.ontology2(results, lfc.cutoff, p.cutoff, correction, direction = 'up', fdr.cut = fdr.cut) %>%
      add_column(Direction = 'up regulated')
    
    go.results.down  = gene.ontology2(results, lfc.cutoff, p.cutoff, correction, direction = 'down', fdr.cut = fdr.cut) %>%
      add_column(Direction = 'down regulated')
    
    go.results <- bind_rows ( go.results.up, go.results.down)
  }
  }
  return(go.results)
}



##########



# eigen plot
# Function to summarize module info for WGCNA
eigenPlot<- function(which.module,low.col,plot.theme = dark_theme_grey()){
  require(tidyverse)
  require(ggdark)
  fill.col<- which.module
  colorh1=moduleColors
  de.p<- t(scale(datExpr[ ,colorh1==which.module]))
  
  de.p %>%
    t() %>%
    as_tibble(rownames = "Genotype") %>%
    mutate(Genotype=sub("(.*-\\d*-\\d*hr).*","\\1", Genotype,perl = T))%>%
    gather(Gene,Expression,-Genotype) %>%
    mutate(Genotype=factor(Genotype)) -> de.p
  
  # set plotting info
  ar<- arrow(angle = 90, length = unit(0.1, "inches"),
             ends = "both" )
  
  H<- ggplot(de.p, aes(x=Gene,y=Genotype,fill=Expression))+
    geom_tile()+
    scale_y_discrete(labels = )+
    scale_fill_gradient2(low = "blue", mid = low.col,high = "red")+
    xlab('Gene Expression')+
    ylab('Genotypes')+
    plot.theme+
    theme(axis.ticks.x=element_blank(),
          axis.text.x=element_text(colour='black'),
          panel.border = element_blank(),
          legend.position = "none")
  
  ## eigen plot
  
  datME[, paste0("ME",which.module)] %>%
    enframe(name=NULL,value="Eigen") %>%
    add_column(Genotype=row.names(datME))%>%
    mutate(Genotype=sub("(.*-\\d*-\\d*hr).*","\\1", Genotype,perl = T))%>%
    mutate(Genotype=factor(Genotype)) -> ME
  
  B<- ggplot(ME,aes(x=Genotype,y=Eigen ,fill=fill.col))+
    geom_bar(stat='summary')+
    stat_summary(geom = "bar", fun.y = mean, position = "dodge",width =.5 ) +
    stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
    scale_fill_manual(values=fill.col)+
    coord_flip() + 
    plot.theme+
    ylab('Eigengene expression')+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          axis.line.x = element_line(color='white',arrow =ar))
  
  grid.arrange(H, B, widths=c(1,.3),
               ncol = 2, nrow = 1)
}


convert.2.fb<- function(x,keytype = 'ENTREZID') {
  AnnotationDbi::select(org.Dm.eg.db,keys = x,columns = 'FLYBASE', keytype = keytype)$FLYBASE
}

convert.2.sym<- function(x,keytype = 'FLYBASE') {
  AnnotationDbi::select(org.Dm.eg.db,keys = x,columns = 'SYMBOL', keytype = keytype)$SYMBOL
}

