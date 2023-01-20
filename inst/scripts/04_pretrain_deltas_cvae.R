reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-4.1-bioc-3.13-reticulate-keras-2.4.0-tensorflow-2.4.0-gpu/")
Sys.setenv("CUDA_VISIBLE_DEVICES" = "0" ) # Define visible  GPU devices
ngpus=length(strsplit(Sys.getenv("CUDA_VISIBLE_DEVICES"),",")[[1]])
reticulate::py_config()
library(keras)
K <- backend() # manual add-on
library(tensorflow)
library(Matrix)
library(SummarizedExperiment)

### Read in summarized experiment with raw counts:
se <-readRDS("mouse_matrix_v212_237Kx20K_se.rds")

#### Library normalize and log transform:
pc <- 4
M <- assay(se)
M <- sweep(M, 2, colSums(M), FUN = "/") * 1e+06
ML <- log2(M+pc)
dimnames(ML) <- dimnames(M)
rm(M)



############################################
## Generate augmented corpus of contrasts ##
############################################

### Create random contrasts within series.
### Acceptance according to distribution of rhos in actual CNT-treatm contrasts
### Highly populated series are subsampled.
### Poorly represented series augmented by adding the inverse Deltas
cor.thr.low <- 0.880   # Lower thr. for pair correlation
cor.thr.upp <- 0.996   # Upper thr. for pair correlation
se$batch_id <- paste0(se$molecule_ch1, se$submission_date, se$contact_name) #Unique batch id

set.seed(1)
smpl.rt <- 0.6 # Sampling  N^(-smpl.rt) fraction of pairs from a series with N concordant pairs.
PAIRS <- list()
N <- list()
bc <- 0

hundredth <- round(length( unique(se$batch_id)  ) / 100)
for (batch in unique(se$batch_id) ){
    bc <- bc +1 
    
    if (bc %% hundredth==0) { cat('\014')
        cat(paste0(round(bc / hundredth), '% completed'))
    }
    
    sel <- se$batch_id== batch 
    n <- sum(sel)
    if (n >1 & n<25) { ## For small series all concordant contrasts are added
        CM <- coop::pcor(ML[,sel])
        CM[lower.tri(CM,diag=TRUE)] <- 0
        concordant <- which(CM > cor.thr.low-0.07/n & CM < cor.thr.upp , arr.ind = TRUE )
        N[[batch]] <- nrow(concordant)
        if(nrow(concordant) >0){
            PAIRS[[batch]] <- cbind( se$geo_accession[sel][concordant[,1]] , se$geo_accession[sel][concordant[,2]] )
        }
        if(nrow(concordant) <16){ # add inverse Delta in poorly represented series
            PAIRS[[batch]] <- rbind(PAIRS[[batch]], PAIRS[[batch]][,2:1] )    
        }
    }
    else if (n > 1 & n>=25){ ## For larger series a fraction of all concordant contrasts are added
        CM <- coop::pcor(ML[,sel])
        CM[lower.tri(CM,diag=TRUE)] <- 0
        concordant <- which(CM > cor.thr.low-0.07/n & CM < cor.thr.upp , arr.ind = TRUE )
        cvals <- CM[concordant]
        if(nrow(concordant) >0){
            nrc <- nrow(concordant)
            
            if (nrc > 100){
                pval <- exp(50*cvals)/5e16 # Pval calculation according to exponential distr. of actual CNT-replicates correlation
                pval <- pval/max(pval)
                smpl.sel <- which(pval > runif(length(pval),0,1)  )
                smpl.sel <- sample(smpl.sel,min(length(smpl.sel),round(100+nrc^smpl.rt )  )  )
                concordant <- concordant[smpl.sel, ]
            }
            
            PAIRS[[batch]] <- cbind( se$geo_accession[sel][concordant[,1]] , se$geo_accession[sel][concordant[,2]] )
            if(nrc <16){ # add inverse Delta in poorly represented series
                PAIRS[[batch]] <- rbind(PAIRS[[batch]], PAIRS[[batch]][,2:1] )    
            }
            
        }
    }
}
sum(unlist(lapply(PAIRS,nrow)))
### ~0.9M Contrasts


### Further filter out contrasts:
### Within every submission batch with more than k pairs check the DLT-DLT correlations of those pairs
### Reduce redundant contrasts by clustering and subsampling
### Remove a fraction of contrasts with low variance (particularly in large series)
k <- 32
PAIRS.filt <- lapply( PAIRS,function(P) {
    thr <- max(0.15, 1.0 - ((nrow(P)^0.2)-2)*0.5)  # Threshold for tree cutting when defining clusters of similar contrasts
    if (nrow(P) > k ) { 
        temp1 <- ML[,P[,1]]
        temp2 <- ML[,P[,2]]
        DLT <- temp1-temp2
        VAR <- apply(DLT,2,var)
        set.seed(1)
        rem <- which(VAR < runif(length(VAR),0.01,  0.02+1e-4 * length(VAR) )) # Remove a fraction of low variance contrasts.
        absCOR <- abs(pcor( DLT ))
        diag(absCOR) <- 0
        
        hcl <- hclust(d = as.dist(1 - absCOR), method = "complete")
        ct <- cutree(hcl, h = 1.0-thr )
        TBL <- table(ct)
        filt.ct <- names(TBL)[TBL > 2]
        if (length(filt.ct >1)) {
            for(c in filt.ct ){
                set.seed(1)
                keep <- sample( which(ct==c)   ,  round(sqrt( TBL[c] ))      ) 
                rem <- c(rem, setdiff(  which(ct==c) , keep ) )       # Remove a fraction of similar (redundant) contrasts
            }
            P <- P[-unique(rem),]
        }
    }
    P   
}
)
sum(unlist(lapply(PAIRS.filt,nrow)))
### ~0.6M Contrasts
P <- unique(do.call(rbind,PAIRS.filt))



###### Cleanup
rm(se)
rm(PAIRS)
rm(PAIRS.filt)
gc()


###### Calculate Deltas and Encode context:
ML <- t(ML) 
temp1 <- ML[P[,1],]
temp2 <- ML[P[,2],]

###### Context input using the average of both conditions:
xCN <- (temp1+temp2 )/2
rm(P)
rm(ML)
gc()

###### Encode context:
encoder <- keras::load_model_hdf5("TRAINED_MODELS/ContextEncoder_ARCHS_v212_64_mouse.hdf5")
latent_output <- predict(encoder, list(gene_input=xCN)) #context
rm(encoder)
gc()
contextL <- ncol(latent_output)
rm(xCN)

####### Delta input 
xD <- temp1-temp2
rm(temp1)
rm(temp2)

####### Splitting in training and validation data:
holdback.samples <- sample(1:nrow(xD),10000) 
train_xD <- xD[-holdback.samples,]
test_xD  <- xD[ holdback.samples,]
rm(xD)
train_xC <- latent_output[-holdback.samples,]
test_xC  <- latent_output[ holdback.samples,]
gc()



##########################################
#### Model definition and compilation ####
##########################################


######## Define the cvae model:
ngenes <- ncol(train_x)
if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
source("models/delta_cVAE.R")


### cVAE loss function:
vae_loss <- function(x, x_decoded_mean){
    reconstruction_loss  <-  loss_mean_squared_error(x, x_decoded_mean)
    kl_loss <- -kl_weight*0.5*K$mean(1 + z_log_varD-log_var_prior - K$square(z_meanD)/var_prior - K$exp(z_log_varD)/var_prior, axis = -1L)  # Delta KL loss
    reconstruction_loss + kl_loss
}

#### Custom correlation function to keep track of. 
cor_metric <- function(y_true, y_pred) {
    y_true_dev <- y_true - k_mean(y_true)
    y_pred_dev <- y_pred - k_mean(y_pred)
    r_num <- k_sum(y_true_dev * y_pred_dev)
    r_den <- k_sqrt(k_sum(k_square(y_true_dev)) * 
                        k_sum(keras::k_square(y_pred_dev)))
    r_num / r_den 
}

#### Compile model
cvae %>% compile(
    loss = vae_loss,
    optimizer = "adam",
    metrics = custom_metric("cor",cor_metric)
)









########################################################################################################  
###########################################     TRAINING      ########################################## 
########################################################################################################  

##### Learning rate scheduler: 
burn_in.nepochs <- 10
burn_in_lr <- 1e-5  
batch_size <- 512


lr_schedule <- function(epoch, current_lr) {
    if (epoch <= burn_in.nepochs ) lr <- burn_in_lr
    else if (epoch < 220) lr <- min (2e-4, burn_in_lr + 5e-5 * ( (epoch-burn_in.nepochs)/10) )#Increase lr linearly 
    else if (epoch > 260 ) {lr <- 1e-6 } #Second LR drop  (Cool down)
    else {lr <- 5e-5} #First LR drop
    return(lr)
}


lr_sch <- callback_learning_rate_scheduler(lr_schedule)

##### Early stopping callback:
early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0,
                                          patience = 250, verbose = 0, mode = "auto",
                                          baseline = NULL, restore_best_weights = TRUE)

##### Checkpoint callback
checkpoint <- callback_model_checkpoint(
    "TRAINED_MODELS/cvae_deltas_ARCHS4_v212_pretr_512_64_mouse.hdf5",
    monitor = "val_loss",
    verbose = 0,
    save_best_only = TRUE,
    save_weights_only = TRUE,
    mode = c( "min"),
    period = NULL,
    save_freq = "epoch"
)


###### Training:
nepochs <- 300
history <- cvae %>% fit(
    x=list(delta_input=train_xD,CONTEXT=train_xC),
    y=train_xD, 
    shuffle = TRUE, 
    epochs = nepochs,
    batch_size = batch_size, 
    validation_data=list(list(test_xD,test_xC),test_xD), 
    callbacks = list( early_stopping, lr_sch, checkpoint )
)

#saveRDS(history,"Trained_models/history_cvae_deltas_ARCHS4_v212_pretr_512_64_mouse.rds")






