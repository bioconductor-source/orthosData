###### Store summarized experiments for different contrast collections:
library(SummarizedExperiment)

se <- readRDS("mouse_matrix_v212_237Kx20K_se.rds")
###### Identify indices of putative controls by searching with a dictionary of 
###### commonly used control-keywords against the submission metadata:
control_keywords <- c("ctr","cntr","contr","cnt","ctl","wt_","wt-","wt ","wt\\.","wild",": wt" ,"no treat","no_treat","no tream",
                      "untr","no[nt]-treat","no[nt]_treat","no[nt] treat","no[nt]tr","empty vect","empty_vect","empty-vect",":GFP",": GFP",
                      "no[nt]targ","no[nt]-targ","no[nt]_targ","no[nt] targ", "mock","scramb","scrm","scrb", "dmso","placebo",":none",
                      ": none","uninf","no[nt] inf","no[nt]_inf","no[nt]-inf","MeOH","shNC","sh-NC","parental","si-NC","si_NC","siNC",
                      "sh_NC","sh\\.NC","shNT","sh_NT","sh-NT","sh\\.NT","siRNA.NC","siRNA-NC","siRNA_NC","vehicle",": dox",
                      "GFP","sgNeg","pretreat","pre-treat","beads only",": healthy",": normal", "buffer alone", ": before", "before .+ treat")

#### Look for the control keywords in the "characteristics_ch1" and "title" fields of the ARCHS4 metadata:
control_idx <- c()
Nmatches <- rep(0,ncol(se))
for (w in control_keywords) {
    for (col in  c("characteristics_ch1","title")){
        idx <- grep(w, colData(se)[,col],ignore.case = TRUE)
        control_idx <- c(control_idx,idx)
        Nmatches[idx] <- Nmatches[idx]+1
    }
}
control_idx <- unique(control_idx)

length(unique(colData(se)$series_id)) #17022
series_with_control <- unique(colData(se)$series_id[control_idx]) #11798 
series_without_control <- setdiff(unique(se$series_id), series_with_control  ) #5224



##### Create deltas for series without identified putative controls by assigning as "control" the series average
date_contact <- paste0( colData(se)$submission_date, colData(se)$contact_name, substr(colData(se)$source_name_ch1,1,30)  )
pc <- 4
C1 <- c()
DELTASnct <- list()
i <- 0
for (ss in unique(date_contact)){
    sel <- which(date_contact==ss)
    if (length(sel) < 3  | sum(sel %in% control_idx)>0 | sum(colData(se)$singlecellprobability[sel] > 0.8) > 10    ) {next}
    i <- i +1
    M <- assay(se)[,sel] 
    M <- sweep(M, 2, colSums(M), FUN = "/") * 1e+06
    mn <- rowMeans(M)
    M <- log2(M+pc)
    mn <- log2(mn+pc)
    COR <- cor(M,mn) 
    names(COR) <- colnames(se)[sel]
    C1 <-  c(C1, COR )
    DELTASnct[[i]] <- M-mn
}
DELTASnct<- do.call(cbind, DELTASnct) #77703
colnames(DELTASnct) <- names(C1)


#### Create deltas for series with identified putative controls and no-controls by pairing them according to max cor:
#### In cases where putative controls have different numbers of matches to the control-word dictionary, assign as control(s)
#### the samples with the maximum number of matches to the controls dictionary.
C2 <- c()
CntNames <- c()
DELTAScnt <-list()
i <- 0
tot_ct <- 0
for (ss in unique(date_contact)){
    sel <- which(date_contact==ss) # Indices of the series
    MaxMatches <- max(Nmatches[sel])
    is.cnt <- Nmatches[sel]==MaxMatches # Assign as controls the series samples with the max number of matches to control keywords
    if (length(sel) < 2  | MaxMatches < 1  | sum(is.cnt) == length(sel) | sum(colData(se)$singlecellprobability[sel] > 0.8) > 10 ) {next}
    tot_ct <- tot_ct + sum(is.cnt)
    i <- i +1 # 8984 series fall in this category
    M <- assay(se)[,sel] 
    M <- sweep(M, 2, colSums(M), FUN = "/") * 1e+06
    M <- log2(M+pc)
    Mctr <- M[, is.cnt,drop=FALSE]
    Mnct <- as.matrix(M[,!is.cnt,drop=FALSE])
    colnames(Mnct) <- se$geo_accession[sel[!is.cnt]]
    COR <-   cor( Mnct, Mctr )
    C2v <-  rowMax(COR) 
    DELTAScnt[[i]] <- Mnct - Mctr[,Rfast::rowMaxs(COR),drop=FALSE]
    CntN  <- se$geo_accession[sel[is.cnt]][Rfast::rowMaxs(COR)]
    names(C2v) <- colnames(Mnct)
    C2 <-  c(C2, C2v ) 
    CntNames <- c(CntNames,CntN)
}
DELTAScnt<- do.call(cbind, DELTAScnt) #71282
colnames(DELTAScnt) <- names(C2)


##### Summarized experiment for all contrasts (with or without an assigned control)
ALL_Deltas <- cbind(DELTASnct,DELTAScnt)
ALL_CORS <- as.data.frame(c(C1,C2))
rownames(ALL_CORS) <- colnames(ALL_Deltas)
HasCNT <- as.data.frame(c( rep(0,length(C1)), rep(1,length(C2))     ))
AllCntNames <-  as.data.frame(c( rep(NA,length(C1)), CntNames    ))
rownames(HasCNT) <- colnames(ALL_Deltas)
DeltaColData <-  cbind( colData(se)[colnames(ALL_Deltas),],ALL_CORS,HasCNT,AllCntNames )
colnames(DeltaColData)[(ncol(DeltaColData)-2):ncol(DeltaColData)] <- c("Cor2CNT","HasAssignedCNT","CNTname")
sed <- SummarizedExperiment(assays=SimpleList(deltas=ALL_Deltas ),
                            rowData=rowData(se),
                            colData=DeltaColData,
                            checkDimnames=TRUE)

sed$aligned_reads <- as.numeric(sed$aligned_reads)
sed$singlecellprobability <- as.numeric(sed$singlecellprobability)
#if (!file.exists("mouse_matrix_v212_149Kx20K_DELTAS_se.rds")) {
#    saveRDS(sed,"mouse_matrix_v212_149Kx20K_DELTAS_se.rds") 
#}


### Filter only contrasts with (reasonably similar) assigned controls
use <- sed$Cor2CNT > 0.85 & sed$Cor2CNT < 0.996 & sed$singlecellprobability < 0.5  & sed$HasAssignedCNT==1 
sedf <- sedf[,use]

### Gather all contexts (for control and perturbations) for the filtered contrasts:
seLCPM_cnt <- seRAW [,sedf$CNTname]
seLCPM <- seRAW [,colnames(sedf)]

###### Context using the average of both conditions:
xC <- assays(seLCPM)[[1]] + assays(seLCPM_cnt)[[1]]
xC <- sweep(xC, 2, colSums(xC), FUN = "/") * 1e+06
xC <- log2(xC+pc)
assays(sedf)[["CONTEXT"]] <- xC

if (!file.exists("mouse_matrix_v212_60Kx20K_DELTASplusCONTEXT_se.rds")) {
    saveRDS(sedf,"mouse_matrix_v212_60Kx20K_DELTASplusCONTEXT_se.rds") 
}





