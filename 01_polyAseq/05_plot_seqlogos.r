library(openxlsx)
library(data.table)
library(ggplot2)
library(ggseqlogo)
library(ggforce) 
library(JASPAR2018)
library(TFBSTools)

fn <- './04_wkd/perm_10K_all.csv'

perm_10k_all_dt <- fread(fn)
perm_10k_all_dt <- perm_10k_all_dt[order(factor)]

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- perm_10k_all_dt$factor

PFMatrixList <- getMatrixSet(JASPAR2018, opts)

PFMatrixList <- PFMatrixList[order(names(PFMatrixList))]

M <- length(PFMatrixList)
N2 <- 10
ncol2 <- 2

pfmatrix <- list()
buffer2 <- 0
pdf_fn <- sprintf('./04_wkd/seqlogos_perm_10k_all.pdf')
pdf(pdf_fn)

for (i in 1:M){
  message(i)
  if (buffer2 >= N2) {
    p <- ggplot() + geom_logo(pfmatrix) + theme_logo() + 
      facet_wrap(~seq_group, ncol=ncol2, scales='free') 
    
    print(p)
    
    buffer2 <- 0
    pfmatrix <- list()
  }
  buffer2 <- buffer2 + 1

  tag <- sprintf('%s(%s)',PFMatrixList[[i]]@name,PFMatrixList[[i]]@ID)

  pfmatrix[[tag]] <- PFMatrixList[[i]]@profileMatrix
  
}

p <- ggplot() + geom_logo(pfmatrix) + theme_logo() + 
  facet_wrap(~seq_group, ncol=ncol2, scales='free') 

print(p)
dev.off()