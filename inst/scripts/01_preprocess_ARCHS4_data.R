###### # Preprocess ARCHS4 gene files to filter data,annotate genes, identify and store contrasts.
library(rhdf5) 
library(Rfast) 
library(SummarizedExperiment)
library(ensembldb)
library(AnnotationHub)


## h5 files downloaded from: https://maayanlab.cloud/archs4/download.html :
## wget https://s3.dev.maayanlab.cloud/archs4/archs4_gene_mouse_v2.1.2.h5 #Mouse
## wget https://s3.dev.maayanlab.cloud/archs4/archs4_gene_human_v2.1.2.h5 #Human
h5_file <-"archs4_gene_mouse_v2.1.2.h5"


# Retrieve information from compressed data
samples <-  h5read(h5_file, "meta/samples/")
samplesDF <- as.data.frame(do.call(cbind,samples))
genes <-  h5read(h5_file, "meta/genes/")
genesDF <-  as.data.frame(do.call(cbind,genes))


##### Filter Samples (select for molecule_type, library_selection, source,
#####                 remove low aln_reads, remove single cell datasete):
molecule.select <- which(samplesDF$molecule_ch1=="polyA RNA" | samplesDF$molecule_ch1=="total RNA" )
library_selection.select <- which(samplesDF$library_selection=="cDNA")
library_source.select <-  which(samplesDF$library_source=="transcriptomic")
aln_reads.select <- which(as.numeric(samplesDF$aligned_reads) > 5e+06)
sc_prob.select <- which(samplesDF$singlecellprobability < 0.5)
select.samples <- Reduce(intersect, list(molecule.select, library_selection.select, library_source.select, aln_reads.select, sc_prob.select)) #236680


##### Retrieve expression data from the hdf5 file:
exp <- h5read(h5_file, "data/expression") 
exp <- exp[select.samples, ]
colnames(exp) <- genesDF$gene_symbol
rownames(exp) <- samplesDF$geo_accession[select.samples]

##### Retrieve Ensdb gene annotation for the ARCHS4 Ensembl version used:
ah <- AnnotationHub()
ENS <- query(ah, c("EnsDb.Mmusculus.v107")) # Mouse
#ENS <- query(ah, c("EnsDb.Hsapiens.v107")) # Human
ensDF <- as.data.frame(genes(ENS[[1]]))
ensDF$gene_name[nchar(ensDF$gene_name) < 2] <- ensDF$gene_id[nchar(ensDF$gene_name) < 2]
ensDF$symbol[nchar(ensDF$symbol) < 2] <- ensDF$gene_id[nchar(ensDF$symbol) < 2]

##### Remove pseudogenes, ribosomal protein genes, genes from non-canonical chromosomes
### Mouse:
pseudog <- grep("pseudo",ensDF$gene_biotype)
noncanon <- unique(c(grep("GL",ensDF$seqnames),grep("JH",ensDF$seqnames),grep("MT",ensDF$seqnames)))
rpgenes <- grep("^M*rp[sl]",ensDF$symbol)
removeg <- unique(c(pseudog, noncanon, rpgenes))

### Human:
#pseudog <- grep("pseudo",ensDF$gene_biotype)
#noncanon <- unique(c(grep("CHR",ensDF$seqnames),grep("GL",ensDF$seqnames),grep("KI",ensDF$seqnames),
#                     grep("LRG",ensDF$seqnames),grep("MT",ensDF$seqnames)))
#rpgenes <- grep("^M*RP[SL]",ensDF$symbol)
#LRGgenes <- grep("^LRG_",ensDF$gene_id)
#removeg <- unique(c(pseudog, noncanon, rpgenes,LRGgenes))

ensDF <- ensDF[-removeg,]
ensDF <- ensDF[!duplicated(ensDF$symbol),]
keep.genes <- genesDF$gene_symbol %in%  ensDF$gene_id | genesDF$gene_symbol %in%  ensDF$symbol 
exp <- exp[,keep.genes]



#### Only keep genes that are among top 10K in at least 2% of the samples OR top 1K in at least 0.2% of the samples:
# Subsample ARCHS4 expression data to speed up calculation:
set.seed(1)
smpl <- sample(1:nrow(exp),50000)
expSmall_t <- t( exp[smpl,] )

rankMat <- nrow(expSmall_t)-Rfast::colRanks(expSmall_t, method="average") #The Rfast descending=TRUE seems buggy
high_expr1 <- which( rowSums(rankMat <= 10000) > ncol(rankMat) * 0.02 )
high_expr2 <- which( rowSums(rankMat <= 1000) > ncol(rankMat) * 0.002 )
high_expr <- unique(c(high_expr1,high_expr2 ))
length(high_expr)
#20339
exp <- exp[,high_expr]


#### Expand ARCHS4 gene annotation using EnsDB data:
ens.idx <- match(colnames(exp),ensDF$gene_id)
ens.idx[is.na(ens.idx)] <- match(colnames(exp)[is.na(ens.idx)],ensDF$gene_name)
RowData <- ensDF[ens.idx,]
rownames(RowData) <- colnames(exp)
col.idx <- match(rownames(exp),samplesDF$geo_accession)
ColData <- samplesDF[col.idx,]
rownames(ColData) <- rownames(exp)

##### Create summarized experiment for filtered expression matrix
se <- SummarizedExperiment(assays=SimpleList(counts=t(exp )),
                           rowData=RowData,
                           colData=ColData,
                           checkDimnames=TRUE)

if (!file.exists("mouse_matrix_v212_237Kx20K_se.rds")) {
    saveRDS(se,"mouse_matrix_v212_237Kx20K_se.rds")
}




