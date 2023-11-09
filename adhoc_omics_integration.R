# This code integrates RNA and ATAC seq data for shARF6 treated cells.
# MIA PaCa-2 cells were treated with shNT or shARF6, and data collected at four time points i.e. days 3, 7, 14 and 21.
# The script is divided into a Main section which has the main function called TF_total_crossref(), that interiorly calls functions that are defined prior to that, namely checkMotif(), atac_open() and atac_closed(). Simply calling the TF_total_crossref() function is sufficient to generate the results presented in the publication.


# Motif checking function -------------------------------------------------

library(tidyverse)

# This data was downloaded from the UCSC web server.
hgTables <- read.table("Data_UCSC/hgTables.txt", header=F) %>% 
  set_names("chr","start","end", "id", "gene")

hgTables$gene <- toupper(hgTables$gene)

# tf=tail(snd$hgnc_symbol, n=1)
# target_list=head(snd$hgnc_symbol, n=-1)

checkMotif <- function(tf,  target_list){
  
  library(MotifDb)
  library(seqLogo)
  library(motifStack)
  library(Biostrings)
  library(GenomicFeatures)
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  
  
  
  # extract MOTIF of TF
  motifs <- query(MotifDb, andStrings=c(tf, "hsapiens"),
                  orStrings=c("jaspar2018", "HOCOMOCOv11"),
                  notStrings = c("::"))
  
  if(length(motifs) == 0){
    
    return(outlist<-list())
    
  } else {
    
    seqLogo(motifs[[1]])
    
    targets <- left_join(as.data.frame(toupper(target_list)) %>% set_names("gene"), hgTables, by="gene") %>% 
      mutate(bed = paste(chr, ":", start, "-", end, sep="")) %>%  na.omit() 
    
    targets = targets[!duplicated(targets$gene),]
    
    modT <- GenomicRanges::makeGRangesFromDataFrame(targets, keep.extra.columns = T)
    promoter.seqs <- getPromoterSeq(modT, Hsapiens, upstream=2000, downstream=300)
    
    pcm.tf <- round(100 * motifs[[1]])
    
    outlist <- lapply(promoter.seqs, matchPWM, pwm=pcm.tf, min.score="70%")
  }
  
  return(outlist)
}

# Returns differential open logFC -----------------------------------
atac_open <- function(gene_list, treat){
  
  alldf <- purrr::map_dfr(.x = list.files("ATAC/Homer") %>% str_subset("UP"),
                          ~ glue::glue("ATAC/Homer/{.x}") %>% read.delim(header = T) %>% dplyr::select(logFC, Gene.Name) %>% 
                            mutate(dataset = .x %>% str_remove("mod_finalresUP.") %>% str_remove(".Homer.txt")) %>% 
                            separate(dataset, into=c("treat","tp")))
  
  open_list <- alldf %>% 
    dplyr::filter(treat=="A6") %>% 
    dplyr::select(-treat) %>% 
    pivot_wider(names_from = tp, values_from = logFC, values_fill=0.000001, values_fn = max) %>% 
    dplyr::select("Gene.Name","3","7","14","21")
  
  
  return(left_join(as.data.frame(toupper(gene_list)) %>% set_names("Gene.Name"), 
                   open_list, 
                   by="Gene.Name") %>% 
           column_to_rownames(var="Gene.Name"))
  
}



# Returns differential closed logFC -----------------------------------
atac_closed <- function(gene_list, treat){
  
  alldf <- purrr::map_dfr(.x = list.files("ATAC/Homer") %>% str_subset("DOWN"),
                          ~ glue::glue("ATAC/Homer/{.x}") %>% read.delim(header = T) %>% dplyr::select(logFC, Gene.Name) %>% 
                            mutate(dataset = .x %>% str_remove("mod_finalresDOWN.") %>% str_remove(".Homer.txt")) %>% 
                            separate(dataset, into=c("treat","tp")))
  
  closed_list <- alldf %>% 
    dplyr::filter(treat=="A6") %>% 
    dplyr::select(-treat) %>% 
    pivot_wider(names_from = tp, values_from = logFC, values_fill=0.00001, values_fn = max) %>% 
    dplyr::select("Gene.Name","3","7","14","21")
  
  return(left_join(as.data.frame(toupper(gene_list)) %>% set_names("Gene.Name"), 
                   closed_list, 
                   by="Gene.Name") %>% 
           column_to_rownames(var="Gene.Name"))
  
}




# Main --------------------------------------------------------------------
TF_total_crossref <- function(treat){
  
  library(tidyverse)
  library(rje)
  library(ensembldb)
  library(EnsDb.Hsapiens.v79)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(reshape)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(stringr)
  library(NbClust)
  library(raster)
  library(dendsort)
  library(viridis)
  
  i = "R6"
  allTg <- data.frame()
  tp<-c("3","7","14","21")
  
  #Transcription factor enrichment was manually performed on the web-based tool for ChEA3 available at https://maayanlab.cloud/chea3/. This data set was downloaded and is available within the RNA/Chea3_R6_hits directory. Here we merge them to create one data frame of TFs and their target genes.
  
  tfdf_scores <- left_join(tidyr::crossing(c(paste("RNA/Chea3_",i,"_hits/Integrated_meanRank_",sep="")),tp,i,"_UP.tsv") %>%  
                             set_names (c("a","b","c","d"))%>%  
                             mutate(names = paste0(a,b,c,d)) %>% 
                             pull(names) %>% 
                             set_names(.,.) %>%  
                             purrr::imap_dfr(.x = ., ~read.delim(.x) %>% 
                                               dplyr::select(TF, Score) %>% 
                                               mutate(file = .y %>% str_remove_all(paste("RNA/Chea3_",i,"_hits/Integrated_meanRank_", sep="")) %>% 
                                                        str_remove_all(".tsv") )) %>% 
                             pivot_wider(names_from =file, values_from = Score),
                           
                           # Finding the overlapping genes
                           tidyr::crossing(c(paste("RNA/Chea3_",i,"_hits/Integrated_meanRank_",sep="")),tp,i,"_UP.tsv") %>%  
                             set_names (c("a","b","c","d"))%>%  
                             mutate(names = paste0(a,b,c,d)) %>% 
                             pull(names) %>% 
                             set_names(.,.) %>%  
                             purrr::imap_dfr(.x = ., ~read.delim(.x) %>% 
                                               dplyr::select(TF, Overlapping_Genes) %>% 
                                               mutate(file = .y %>% str_remove_all(paste("RNA/Chea3_",i,"_hits/Integrated_meanRank_", sep="")) %>% 
                                                        str_remove_all(".tsv") )),
                           # %>% 
                           #   pivot_wider(names_from =file), 
                           by="TF")
  
  #Top TFs are selected based on the ChEA3 score. Here it has been filtered for <100
  
  tfdf_scores <- cbind(tfdf_scores, min_score = rowMins(tfdf_scores[,2:5])) %>% 
    dplyr::filter(min_score < 100)
  
  col.order <- c("TF",
                 paste("3",i,"_UP", sep=""),
                 paste("7",i,"_UP", sep=""),
                 paste("14",i,"_UP", sep=""),
                 paste("21",i,"_UP", sep=""),
                 "Overlapping_Genes")
  
  tfdf_scores <- tfdf_scores[,col.order]
  
  
  # Compiling the target genes of the TFs.
  
  targets <- list()
  
  for(m in 1:length(unique(tfdf_scores$TF))*4){
    targets[[m]] <- as.data.frame(str_split(paste(tfdf_scores$Overlapping_Genes[m-3],
                                                  tfdf_scores$Overlapping_Genes[m-2],
                                                  tfdf_scores$Overlapping_Genes[m-1],
                                                  tfdf_scores$Overlapping_Genes[m], sep=","), ","))  %>%  
      
      set_names("target_genes") %>%  
      distinct()
    
  }
  
  targets <- targets[lengths(targets) != 0]
  names(targets) <- tfdf_scores$TF %>% unique()
  
  tfdf_scores <- tfdf_scores %>% 
    dplyr::select(1:5) %>% 
    distinct() %>% 
    column_to_rownames(var="TF")
  
  
  # Normalised counts data was used from RNAseq
  
  normcount_data=read.delim("RNA/normalized_counts_log2_star.txt", header=T)
  ensembl_genenames=str_replace(rownames(normcount_data),pattern=".[0-9]+$",replacement="")
  hgnc_genenames <- ensembldb::select(EnsDb.Hsapiens.v79, keys=ensembl_genenames, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
  normcount_data$hgnc_symbol=cbind(ensembl_genenames,hgnc_genenames[match(ensembl_genenames,hgnc_genenames$GENEID),1])[,2]  
  
  pb <- winProgressBar(title = "progress bar", min = 0,max = length(names(targets)), width = 300)
  
  # Loop through all tfs
  for(k in 1:length(names(targets))){
    print(k) 
    
    # Create a sub df of target genes and respective TFs
    gene_list <- rbind(targets[[k]], names(targets)[k]) %>%  set_names("hgnc_symbol")
    snd <- suppressMessages(left_join(gene_list, normcount_data)) %>% na.omit()
    
    snd <- snd[!duplicated(snd$hgnc_symbol, fromLast=T),]
    
    # Checkpoint to exit if TF is not found in the normcount_data file, i.e. if the values are NA.
    if(names(targets)[k] %in% snd$hgnc_symbol == FALSE){
      print(paste(names(targets)[k], " was not found in the normacount_data", sep=""))
      next
    }  
    
    
    # Check for presence or absence of TF motifs on promoters of target genes ---------------------------------------
    
    checkM <- checkMotif(tail(snd$hgnc_symbol, n=1),  head(snd$hgnc_symbol, n=-1))
    
    check <- checkM %>% purrr::map_lgl(.x = ., ~(.x %>% length()) > 0) %>% 
      as.data.frame() %>% 
      set_names('motif_check')%>% 
      rownames_to_column(var = 'index') %>% 
      mutate(hgnc_symbol = snd$hgnc_symbol[as.integer(index)]) %>% 
      dplyr::select(hgnc_symbol, motif_check) %>%
      na.omit()
    
    check <- left_join(gene_list, check, by="hgnc_symbol") %>% 
      set_names("gene_name", "checkOut")
    
    # check if checkMotif() failed i.e. returned ONLY NAs. If yes, then exit
    if(is.element(TRUE, check$checkOut)==FALSE && is.element(FALSE, check$checkOut)==FALSE){ 
      next
    } 
    
    
    #Creating the normalized list
    mean_exp <- list( 
      
      datRN=(data.frame(snd[,2:13])) %>% 
        rowwise() %>% 
        mutate(T3 = mean(c(RN0301, RN0302, RN0303)),
               T7 = mean(c(RN0701, RN0702, RN0703)),
               T14 = mean(c(RN1401, RN1402, RN1403)),
               T21 = mean(c(RN2101, RN2102, RN2103))),
      datRK=(data.frame(snd[,14:19])) %>%
        rowwise() %>%
        mutate(T3 = mean(c(RK0301, RK0302, RK0303)),
               T7 = mean(c(RK0701, RK0702, RK0703))),
      datR6=(data.frame(snd[,26:37])) %>% 
        rowwise() %>% 
        mutate(T3 = mean(c(R60301, R60302, R60303)),
               T7 = mean(c(R60701, R60702, R60703)),
               T14 = mean(c(R61401, R61402, R61403)),
               T21 = mean(c(R62101, R62102, R62103))),
      datR1=(data.frame(snd[,38:49])) %>%
        rowwise() %>%
        mutate(T3 = mean(c(R10301, R10302, R10303)),
               T7 = mean(c(R10701, R10702, R10703)),
               T14 = mean(c(R11401, R11402, R11403)),
               T21 = mean(c(R12101, R12102, R12103)))
    )
    
    
    plot_p <- data.frame()
    plot_p_tf <- data.frame()
    anno_tf_enriched <- data.frame()
    
    plot_p <- cbind(as.data.frame((mean_exp$datR6-mean_exp$datRN)[,13:16]),
                    gene_name=snd$hgnc_symbol) %>% 
      left_join(check, by="gene_name") %>% 
      distinct() %>% 
      column_to_rownames(var="gene_name")
    
    plot_p_tf <- tail(plot_p, n=1)[,1:4] # TF expression levels for top annotation
    
    plot_p_motifFound <- head(plot_p,-1) %>% 
      dplyr::filter(checkOut==TRUE)
    
    
    anno_tf_enriched <- tfdf_scores[names(targets)[k],] # TF motif enrichment scores for top annotation
    colnames(anno_tf_enriched) <- colnames(plot_p)[1:4]
    
    
    if(nrow(plot_p_motifFound) < 10 || nrow(plot_p_tf) == 0 || nrow(anno_tf_enriched) == 0){
      print(paste(names(targets)[k],"was eliminated after filtering", sep=" "))
      next
    }
    
    # filter genes for those which directionally agree with the TF
    trend <- plot_p_motifFound[1:4] / plot_p_tf[col(plot_p_motifFound[1:4])] 
    trend <- rownames_to_column(trend, var = 'gene') %>% 
      rowwise() %>% 
      mutate(Min = min(T3,T7,T14,T21)) %>% 
      dplyr::filter(Min > 0) %>% 
      dplyr::select(gene)
    
    plot_p_motifFound <- left_join(trend, 
                                   plot_p_motifFound %>% 
                                     rownames_to_column(var="gene"),
                                   by="gene") %>% 
      column_to_rownames("gene")
    
    if(nrow(plot_p_motifFound) < 10 || nrow(plot_p_tf) == 0 || nrow(anno_tf_enriched) == 0){
      print(paste(names(targets)[k],"was eliminated after filtering", sep=" "))
      next
    }
    
    # Prepare atac logFC peak open/close ----------------------------------------------
    checkA_open <- atac_open(rownames(plot_p_motifFound),i)
    checkA_open[is.na(checkA_open)] <- 0.00001
    
    checkA_closed <- atac_closed(rownames(plot_p_motifFound),i)
    checkA_closed[is.na(checkA_closed)] <- 0.00001
    
    checkA <- left_join(as.data.frame(rownames(plot_p_motifFound)) %>% set_names("Gene.Name"), 
                        
                        rbind(checkA_closed%>% rownames_to_column(var = "Gene.Name"),
                              checkA_open%>% rownames_to_column(var = "Gene.Name")) %>% 
                          distinct() %>%
                          group_by(Gene.Name) %>%
                          summarise(T3 = `3`[which.max(abs(`3`))],
                                    T7 = `7`[which.max(abs(`7`))],
                                    T14 = `14`[which.max(abs(`14`))],
                                    T21 = `21`[which.max(abs(`21`))]),
                        # purrr::map(.x = .,~.x %>% dplyr::select(where(is.numeric)) %>% colMeans(na.rm = T)),
                        
                        by = "Gene.Name") %>% 
      column_to_rownames("Gene.Name")
    
    
    hc <- hclust(dist(plot_p_motifFound[,1:4]), method="complete")
    
    
    nb <- NbClust(data = plot_p_motifFound[,1:4], distance = "euclidean", min.nc = 1, max.nc = 9, 
                  method = "complete", index = "kl", alphaBeale = 0.1)
    
    # sp <- sp_corr(x=plot_p_motifFound[,1:4], y=checkA)
    sp_corr <- data.frame()
    for(s in 1:length(rownames(plot_p_motifFound))){
      sp_corr <- rbind(sp_corr, cor(as.numeric(plot_p_motifFound[s,1:4]),
                                    as.numeric(checkA[s,1:4]), method="spearman"))
    }
    sp_corr <- sp_corr %>% set_names('sp')
    
    plot_df <- cbind(plot_p_motifFound, checkA, sp_corr, nb$Best.partition) %>% 
      set_names("rna_t3","rna_t7","rna_t14","rna_t21", "motif", 
                "atac_t3","atac_t7","atac_t14","atac_t21",
                "sp", "cluster")
    
    cluster_filter <- plot_df %>% 
      group_by(cluster) %>% 
      summarise(n=median(sp, na.rm=T)) %>% 
      dplyr::filter(n>0) %>% 
      dplyr::select(cluster)
    
    plot_df <- plot_df %>% 
      dplyr::filter(cluster %in% cluster_filter$cluster)
    
    
    if(nrow(plot_df) < 10 || nrow(plot_p_tf) == 0 || nrow(anno_tf_enriched) == 0){
      print(paste(names(targets)[k],"was eliminated after sp_corr", sep=" "))
      next
    }
    
    # heatmaps ----------------------------------------------------------------
    # heatmap colors
    # col_rnorm = colorRamp2(c(-2, 0, 2), c("#0077b6", "#FAF9F6", "#d00000"))
    col_rnorm = colorRamp2(c(-2, 0, 2), c(viridis(100)[1],viridis(100)[50],viridis(100)[100]))
    col_MO = colorRamp2(c(400, 0), c( "#FAF9F6", "#1fc600"))
    col_atac = colorRamp2(c(-2, 0, 2), c("#1fc600", "#FAF9F6", "#e36414"))
    
    
    col_rnorm = colorRamp2(c(-2, 0, 2), c(viridis(100)[1],viridis(100)[50],viridis(100)[100]))
    
    
    
    # heatmap 1 annotations
    ha_top <- HeatmapAnnotation(TF = plot_p_tf %>% unlist(),
                                MO = anno_tf_enriched %>% unlist(),
                                col = list(TF = col_rnorm, MO = col_MO),
                                show_legend = FALSE)
    ha_top@anno_list$TF@label <- names(targets)[k]
    ha_top@anno_list$MO@label <- "TF enrichment"
    
    # right annotation is a check for MOTIF presence in teh promoter region of each gene
    ha_right <- rowAnnotation(motifCheck = plot_df$motif,
                              col = list(motifCheck = c('TRUE'="purple",'FALSE'="black")),
                              show_legend = FALSE)
    
    # heatmap 1
    if(length(unique(plot_df$cluster)) > 1){
      
      set.seed(171)
      ht1 <- Heatmap(plot_df[,1:4],
                     cluster_columns = F,
                     row_split = length(unique(plot_df$cluster)),
                     cluster_rows = dendsort(hclust(dist(plot_df[,1:4]), method="complete")),
                     show_row_names = F,
                     col = col_rnorm,
                     top_annotation = ha_top,
                     right_annotation = ha_right,
                     row_gap = unit(2, "mm"),
                     width = unit(7, "cm"),
                     show_heatmap_legend = T)
    } else {
      
      set.seed(171)
      ht1 <- Heatmap(plot_df[,1:4],
                     cluster_columns = F,
                     cluster_rows = dendsort(hclust(dist(plot_df[,1:4]), method="complete")),
                     show_row_names = F,
                     col = col_rnorm,
                     top_annotation = ha_top,
                     right_annotation = ha_right,
                     row_gap = unit(2, "mm"),
                     width = unit(7, "cm"),
                     show_heatmap_legend = T)
    }
    
    
    # heatmap 2 annotations
    panel_fun = function(index, nm) {
      pushViewport(viewport(xscale = c(-1,1), yscale = c(0, 2)))
      grid.rect()
      grid.xaxis(gp = gpar(fontsize = 8))
      grid.boxplot(plot_df$sp[index] %>% na.omit(), pos = 1, direction = "horizontal")
      popViewport()
    }
    
    anno = anno_link(align_to = plot_df$cluster, which = "row", panel_fun = panel_fun,
                     size = unit(2, "cm"), gap = unit(1, "cm"), width = unit(4, "cm"))
    set.seed(171)
    ht2 <-  Heatmap(plot_df[,6:9],
                    cluster_columns = F,
                    right_annotation = rowAnnotation(foo = anno),
                    show_row_names = F,
                    row_gap = unit(2, "mm"),
                    width = unit(4, "cm"),
                    show_heatmap_legend = T,
                    col = col_rnorm)
    
    
    
    
    
    # Print heatmaps ----------------------------------------------     
    png(file=paste("Output/",i,"/",names(targets)[k], ".png", sep = ""),
        width = 720, height = 720,
        pointsize = 12)
    
    set.seed(171)
    draw(ht1+ht2) #print both heatmaps together
    
    dev.off()
    
    # Target gene - tf data frame --------------------------------
    allTg <- rbind(allTg,
                   as.data.frame(rownames(plot_df)) %>% 
                     mutate(tf=names(targets)[k]) %>% 
                     set_names("tg", "tf")
    )
    
    # Target gene - tf data frame --------------------------------
    allTg <- rbind(allTg,
                   as.data.frame(rownames(plot_df)) %>% 
                     mutate(tf=names(targets)[k]) %>% 
                     set_names("tg", "tf")
    )
    
    # Update perecntage completion
    percent = ceiling((k/length(names(targets)) * 100))
    setWinProgressBar(pb, k, title=paste(i, ": ", percent,"% done", sep=""))
    
  }
  
  write.csv(allTg, file=paste("Output/",i,"/TG.csv", sep=""),
            col.names=T,
            quote=F,
            row.names = F)
  
  close(pb)
  
}


# Reactome Analysis -------------------------------------------------------

nt_normal <- read.table("RNA/normalized_counts_relativeNT.txt") %>% 
  rownames_to_column(var="elim") %>% 
  dplyr::select(-elim)

# Reading in genes which are the target genes of important TFs
allgenesr6 <- read.csv("Output/R6/TG.csv", header=T) %>% 
  dplyr::select(tg) %>% set_names("symbol") %>% 
  left_join(nt_normal, by="symbol") %>% 
  dplyr::select(symbol,R603,R607,R614,R621) %>% 
  distinct()

up_shortlist <- data.frame()
adup_shortlist <- data.frame()

for(i in 1:length(allgenesr6$symbol)){
  
  if(allgenesr6[i,2]>=0 && allgenesr6[i,3]>=0 && allgenesr6[i,4]>=0 && allgenesr6[i,5] >=0){
    
    up_shortlist <- rbind(up_shortlist, allgenesr6[i,])
    
  } else if(allgenesr6[i,2]<0 && (allgenesr6[i,4]>0 & allgenesr6[i,5]>0)){
    
    adup_shortlist <- rbind(adup_shortlist, allgenesr6[i,])
    
  } 
}

x <- up_shortlist %>% 
  pivot_longer(cols = 2:5, names_to='tp', values_to = 'log2fc') %>% 
  mutate(tp = factor(tp, levels=c("R603", "R607", "R614", "R621")))

ggup <- ggplot(data=x, aes(x=tp, y=log2fc, group=symbol))+
  # geom_line(color="#398AD7", aes(alpha=0.2))+
  geom_point()+
  ylim(-5,5)+
  ggtitle(label="Upregulation")+
  theme_bw()+ theme(legend.position = "none")+ 
  stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1))

y <- adup_shortlist %>% 
  pivot_longer(cols = 2:5, names_to='tp', values_to = 'log2fc') %>% 
  mutate(tp = factor(tp, levels=c("R603", "R607", "R614", "R621")))

ggadup <- ggplot(data=y, aes(x=tp, y=log2fc, group=symbol))+
  # geom_line(color="#398AD7", aes(alpha=0.2))+
  geom_point()+
  ylim(-5,5)+
  ggtitle(label="Adaptive upregulation")+
  theme_bw()+ theme(legend.position = "none")+ 
  stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1))

grid.arrange(ggup, ggadup, ncol=2)

pdf("Output/expression_based_clusters.pdf", width=5, height=2.5)
grid.arrange(ggup, ggadup, ncol=2)
dev.off()

#### Reactome
# The [ReactomePA](https://bioconductor.org/packages/release/bioc/html/ReactomePA.html) package provides functions for pathway analysis based on REACTOME pathway database. It implements enrichment analysis using a main list of genes and a background list. 

r6_up <- up_shortlist %>% dplyr::select(symbol) %>% unique()
r6_adup <- adup_shortlist %>% dplyr::select(symbol) %>% unique()

# function to remove NAs from list
na.omit.list <- function(y) { return(y[!vapply(y, function(x) all(is.na(x)), logical(1))]) }

# entrez -- genes
r6_up_Entrez <- mapIds(org.Hs.eg.db, r6_up$symbol, 'ENTREZID', 'SYMBOL') %>% 
  na.omit.list() %>%
  unlist() %>%
  as.data.frame(row.names = names(.)) %>% 
  distinct() %>% 
  set_names("entrez") 

r6_adup_Entrez <- mapIds(org.Hs.eg.db, r6_adup$symbol, 'ENTREZID', 'SYMBOL') %>% 
  na.omit.list() %>%
  unlist() %>%
  as.data.frame(row.names = names(.)) %>% 
  distinct() %>% 
  set_names("entrez") 


# From the Reactome analysis, we were only able to find enriched terms for either *upregulated* or *adaptively upregulated* groups.

r6_reactome_up <- enrichPathway(
  r6_up_Entrez$entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  # universe = bg_Entrez$entrez,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
r6_reactome_up@result$Description[14] = "TLR7/8 --> TRAF6 --> NFkB/MAPK"

r6_reactome_adup <- enrichPathway(
  r6_adup_Entrez$entrez,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  # universe = bg_Entrez$entrez,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE
)
r6_reactome_adup@result$Description[13] = "LLO transfer to nacent protein"

r6_reactome_plot_up <- data.frame( Description = r6_reactome_up$Description, 
                              pvalue = r6_reactome_up$p.adjust, 
                              Count = r6_reactome_up$Count) %>% 
  head(20) %>% 
  ggplot(aes(x=Count, y=reorder(Description,Count), fill=pvalue))+
  geom_col()+
  scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.grid = element_blank())

r6_reactome_plot_adup <- data.frame( Description = r6_reactome_adup$Description, 
                                     pvalue = r6_reactome_adup$p.adjust, 
                                     Count = r6_reactome_adup$Count) %>% 
  head(20) %>% 
  ggplot(aes(x=Count, y=reorder(Description,Count), fill=pvalue))+
  geom_col()+
  scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.grid = element_blank())

title1 = text_grob("Constitutive upregulation", size = 15, face = "bold") 
title2 = text_grob("Adaptive upregulation", size = 15, face = "bold") 

pdf("Output/ReactomePA.pdf", width=13, height=4)
grid.arrange(arrangeGrob(r6_reactome_plot_up, top = title1),
             arrangeGrob(r6_reactome_plot_adup, top = title2),
             layout_matrix = rbind(c(1,2)))
dev.off()