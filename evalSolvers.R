library(TReNA)
library(plyr)
library(limma)


assess_methodsAgainstDistributions <- function(mtx.sub, target.gene, tfs)
{    
  #fivenum(mtx.sub)# 0.000000    1.753137   12.346965   43.247467 1027.765854
  
  # Transform with log2 
  mtx.tmp <- mtx.sub - min(mtx.sub) + 0.001
  mtx.log2 <- log2(mtx.tmp)
  #fivenum(mtx.log2)  # [1] -9.9657843  0.8107618  3.6262014  5.4345771 10.005297
  
  # Transform sub with asinh
  mtx.asinh <- asinh(mtx.sub)
  #fivenum(mtx.asinh)  # [1] 0.000000 1.327453 3.208193 4.460219 7.628290
  
  # Transform via VOOM transformation
  mtx.voom <- voom(mtx.sub)$E
  #fivenum(mtx.voom) 
  
  print("--- Testing Spearman")
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="spearman", quiet=FALSE)
  spearman1 <- solve(trena, target.gene, tfs)
  spearman1 <- data.frame(gene = rownames(spearman1),
                          spearman.as.is = spearman1$coefficient)
  
  
  trena <- TReNA(mtx.assay=mtx.log2, solver="spearman", quiet=FALSE)
  spearman2 <- solve(trena, target.gene, tfs)
  spearman2 <- data.frame(spearman.log2 = spearman2$coefficient,
                          gene = rownames(spearman2))
  
  trena <- TReNA(mtx.assay=mtx.asinh, solver="spearman", quiet=FALSE)
  spearman3 <- solve(trena, target.gene, tfs)
  spearman3 <- data.frame(spearman.asinh = spearman3$coefficient,
                          gene = rownames(spearman3))
  
  trena <- TReNA(mtx.assay=mtx.voom, solver="spearman", quiet=FALSE)
  spearman4 <- solve(trena, target.gene, tfs)
  spearman4 <- data.frame(spearman.voom = spearman4$coefficient,
                          gene = rownames(spearman4))
  
  
  spearman1$gene <- as.character(spearman1$gene)
  spearman2$gene <- as.character(spearman2$gene)
  spearman3$gene <- as.character(spearman3$gene)
  spearman4$gene <- as.character(spearman4$gene)
  
  # Grab the top 10 genes from each
  spearman1.top10 <- spearman1$gene[order(abs(spearman1$spearman.as.is), decreasing=TRUE)][1:10]
  spearman2.top10 <- spearman2$gene[order(abs(spearman2$spearman.log2), decreasing=TRUE)][1:10]
  spearman3.top10 <- spearman3$gene[order(abs(spearman3$spearman.asinh), decreasing=TRUE)][1:10]
  spearman4.top10 <- spearman4$gene[order(abs(spearman4$spearman.voom), decreasing = TRUE)][1:10]
  
  print("--- Testing Ridge")
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="ridge", quiet=FALSE)
  result.ridge1 <- solve(trena, target.gene, tfs)
  ridge1 <- data.frame(ridge.as.is = result.ridge1$beta,
                       gene = rownames(result.ridge1))
  
  trena <- TReNA(mtx.assay=mtx.log2, solver="ridge", quiet=FALSE)
  result.ridge2 <- solve(trena, target.gene, tfs)
  ridge2 <- data.frame(ridge.log2 = result.ridge2$beta,
                       gene = rownames(result.ridge2))
  
  trena <- TReNA(mtx.assay=mtx.asinh, solver="ridge", quiet=FALSE)
  result.ridge3 <- solve(trena, target.gene, tfs)
  ridge3 <- data.frame(ridge.asinh = result.ridge3$beta,
                       gene = rownames(result.ridge3))
  
  trena <- TReNA(mtx.assay=mtx.voom, solver="ridge", quiet=FALSE)
  result.ridge4 <- solve(trena, target.gene, tfs)
  ridge4 <- data.frame(ridge.voom = result.ridge4$beta,
                       gene = rownames(result.ridge4))
  
  ridge1$gene <- as.character(ridge1$gene)
  ridge2$gene <- as.character(ridge2$gene)
  ridge3$gene <- as.character(ridge3$gene)
  ridge4$gene <- as.character(ridge4$gene)
  
  # Grab the top 10 genes from each
  ridge1.top10 <- ridge1$gene[order(abs(ridge1$ridge.as.is), decreasing=TRUE)][1:10]
  ridge2.top10 <- ridge2$gene[order(abs(ridge2$ridge.log2), decreasing=TRUE)][1:10]
  ridge3.top10 <- ridge3$gene[order(abs(ridge3$ridge.asinh), decreasing=TRUE)][1:10]
  ridge4.top10 <- ridge4$gene[order(abs(ridge4$ridge.voom), decreasing = TRUE)][1:10]
  
  print("--- Testing Lassopv")
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="lassopv", quiet=FALSE)
  lpv.result <- solve(trena, target.gene, tfs)
  lpv1 <- data.frame(lpv.as.is = -log10(lpv.result$p.values),
                     gene = rownames(lpv.result))
  
  trena <- TReNA(mtx.assay=mtx.log2, solver="lassopv", quiet=FALSE)
  lpv.result.2 <- solve(trena, target.gene, tfs)
  lpv2 <- data.frame(lpv.log2 = -log10(lpv.result.2$p.values),
                     gene = rownames(lpv.result.2))
  
  trena <- TReNA(mtx.assay=mtx.asinh, solver="lassopv", quiet=FALSE)
  lpv.result.3 <- solve(trena, target.gene, tfs)
  lpv3 <- data.frame(lpv.asinh = -log10(lpv.result.3$p.values),
                     gene = rownames(lpv.result.3))
  
  trena <- TReNA(mtx.assay=mtx.voom, solver="lassopv", quiet=FALSE)
  lpv.result.4 <- solve(trena, target.gene, tfs)
  lpv4 <- data.frame(lpv.voom = -log10(lpv.result.4$p.values),
                     gene = rownames(lpv.result.4))
  
  lpv1$gene <- as.character(lpv1$gene)
  lpv2$gene <- as.character(lpv2$gene)
  lpv3$gene <- as.character(lpv3$gene)
  lpv4$gene <- as.character(lpv4$gene)
  
  # Grab the top 10 genes from each
  lpv1.top10 <- lpv1$gene[order(abs(lpv1$lpv.as.is), decreasing=TRUE)][1:10]
  lpv2.top10 <- lpv2$gene[order(abs(lpv2$lpv.log2), decreasing=TRUE)][1:10]
  lpv3.top10 <- lpv3$gene[order(abs(lpv3$lpv.asinh), decreasing=TRUE)][1:10]    
  lpv4.top10 <- lpv4$gene[order(abs(lpv4$lpv.voom), decreasing=TRUE)][1:10]
  
  # Use ensemble!
  print("--- Testing Ensemble")
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="ensemble", quiet=FALSE)
  ensemble1 <- solve(trena, target.gene, tfs,
                     extraArgs = list(solver.list = c("lasso",
                                                      "ridge",
                                                      "randomForest",
                                                      "sqrtlasso",
                                                      "lassopv",
                                                      "pearson",
                                                      "spearman")))
  ensemble1 <- data.frame(gene = ensemble1$gene,
                          ensemble.as.is = ensemble1$pcaMax)
  
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="ensemble", quiet=FALSE)
  ensemble2 <- solve(trena, target.gene, tfs,
                     extraArgs = list(solver.list = c("lasso",
                                                      "ridge",
                                                      "randomForest",
                                                      "sqrtlasso",
                                                      "lassopv",
                                                      "pearson",
                                                      "spearman")))
  ensemble2 <- data.frame(ensemble.log2 = ensemble2$pcaMax,
                          gene = ensemble2$gene)
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="ensemble", quiet=FALSE)
  ensemble3 <- solve(trena, target.gene, tfs,
                     extraArgs = list(solver.list = c("lasso",
                                                      "ridge",
                                                      "randomForest",
                                                      "sqrtlasso",
                                                      "lassopv",
                                                      "pearson",
                                                      "spearman")))
  ensemble3 <- data.frame(ensemble.asinh = ensemble3$pcaMax,
                          gene = ensemble3$gene)
  
  trena <- TReNA(mtx.assay=mtx.sub, solver="ensemble", quiet=FALSE)
  ensemble4 <- solve(trena, target.gene, tfs,
                     extraArgs = list(solver.list = c("lasso",
                                                      "ridge",
                                                      "randomForest",
                                                      "sqrtlasso",
                                                      "lassopv",
                                                      "pearson",
                                                      "spearman")))
  ensemble4 <- data.frame(ensemble.voom = ensemble4$pcaMax,
                          gene = ensemble4$gene)
  
  
  ensemble1$gene <- as.character(ensemble1$gene)
  ensemble2$gene <- as.character(ensemble2$gene)
  ensemble3$gene <- as.character(ensemble3$gene)
  ensemble4$gene <- as.character(ensemble4$gene)
  
  # Grab the top 10 genes from each
  ensemble1.top10 <- ensemble1$gene[order(abs(ensemble1$ensemble.as.is), decreasing=TRUE)][1:10]
  ensemble2.top10 <- ensemble2$gene[order(abs(ensemble2$ensemble.log2), decreasing=TRUE)][1:10]
  ensemble3.top10 <- ensemble3$gene[order(abs(ensemble3$ensemble.asinh), decreasing=TRUE)][1:10]
  ensemble4.top10 <- ensemble4$gene[order(abs(ensemble4$ensemble.voom), decreasing = TRUE)][1:10]
  
  
  # Take the union of all the genes
  all.genes <- unique(c(spearman1.top10,
                        spearman2.top10,
                        spearman3.top10,
                        spearman4.top10,
                        ridge1.top10,
                        ridge2.top10,
                        ridge3.top10,
                        ridge4.top10,
                        lpv1.top10,
                        lpv2.top10,
                        lpv3.top10,
                        lpv4.top10,
                        ensemble1.top10,
                        ensemble2.top10,
                        ensemble3.top10,
                        ensemble4.top10)) # Returns merged genes (37 total)
  
  
  
  # Pull out the specified genes and assemble a table
  spearman1.sub <- subset(spearman1, spearman1$gene %in% all.genes)
  spearman2.sub <- subset(spearman2, spearman2$gene %in% all.genes)
  spearman3.sub <- subset(spearman3, spearman3$gene %in% all.genes)
  spearman4.sub <- subset(spearman4, spearman4$gene %in% all.genes)
  ridge1.sub <- subset(ridge1, ridge1$gene %in% all.genes)
  ridge2.sub <- subset(ridge2, ridge2$gene %in% all.genes)
  ridge3.sub <- subset(ridge3, ridge3$gene %in% all.genes)
  ridge4.sub <- subset(ridge4, ridge4$gene %in% all.genes)
  lpv1.sub <- subset(lpv1, lpv1$gene %in% all.genes)
  lpv2.sub <- subset(lpv2, lpv2$gene %in% all.genes)
  lpv3.sub <- subset(lpv3, lpv3$gene %in% all.genes)
  lpv4.sub <- subset(lpv4, lpv4$gene %in% all.genes)
  ensemble1.sub <- subset(ensemble1, ensemble1$gene %in% all.genes)
  ensemble2.sub <- subset(ensemble2, ensemble2$gene %in% all.genes)
  ensemble3.sub <- subset(ensemble3, ensemble3$gene %in% all.genes)
  ensemble4.sub <- subset(ensemble4, ensemble4$gene %in% all.genes)
  
  # Join it all in a table (requires plyr package)
  tbl.all <- join_all(list(spearman1.sub,
                           spearman2.sub,
                           spearman3.sub,
                           spearman4.sub,
                           ridge1.sub,
                           ridge2.sub,
                           ridge3.sub,
                           ridge4.sub,
                           lpv1.sub,
                           lpv2.sub,
                           lpv3.sub,
                           lpv4.sub,
                           ensemble1.sub,
                           ensemble2.sub,
                           ensemble3.sub,
                           ensemble4.sub),
                      by = 'gene', type = 'full')
  
  # Replace NA's in any columns with 0s
  tbl.all[is.na(tbl.all)] <- 0
  
  # Add the gene correlation column
  tbl.all$gene.cor <- NA
  for(gene in tbl.all$gene){
    tbl.all$gene.cor[[which(tbl.all$gene == gene)]] <-
      mean(result.ridge1$gene.cor[[which(rownames(result.ridge1) == gene)]],
           lpv.result$gene.cor[[which(rownames(lpv.result) == gene)]])
  }
  
  # Order the rows and return it 
  tbl.all <- tbl.all[order(abs(tbl.all$gene.cor), decreasing = TRUE),]
  invisible(tbl.all)
  
} #assess_methodsAgainstDistributions
#----------------------------------------------------------------------------------------------------
assess_ampAD154AllSolversAndDistributions <- function(){
  
  # Load the matrix and create the transformed versions
  load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
  target.gene <- "MEF2C"
  tfs <- setdiff(rownames(mtx.sub), "MEF2C")
  
  tbl.all <- assess_methodsAgainstDistributions(mtx.sub,target.gene,tfs)
  
  
} #assess_ampAD154AllSolversAndDistributions
#----------------------------------------------------------------------------------------------------