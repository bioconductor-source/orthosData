reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-4.1-bioc-3.13-reticulate-keras-2.4.0-tensorflow-2.4.0-gpu/")
Sys.setenv("CUDA_VISIBLE_DEVICES" = "0" ) # Define visible  GPU devices

ngpus=length(strsplit(Sys.getenv("CUDA_VISIBLE_DEVICES"),",")[[1]])
reticulate::py_config()
library(keras)
K <- backend() # manual add-on
library(tensorflow)
library(Matrix)
library(SummarizedExperiment)


## Load in summarized experiment of contrasts with valid assigned controls (includes DELTAS + CONTEXT assays)
se <- readRDS("mouse_matrix_v212_60Kx20K_DELTASplusCONTEXT_se.rds")    
ngenes <- nrow(se)

###### Context input using the average of both conditions:
xC <- assays(se)[["CONTEXT"]]
xC <- t(xC)                                 #Need to transpose before passing to the model
####### Delta input 
xD <- assays(se)[["deltas"]]
xD <- t(xD)                                 #Need to transpose before passing to the model

###### Encode context:
encoder <- keras::load_model_hdf5("TRAINED_MODELS/ContextEncoder_ARCHS_v212_64_mouse.hdf5")
latent_output <- predict(encoder, list(gene_input=xC)) #context
rm(encoder)
gc()
contextL <- ncol(latent_output)


####### Splitting in training and validation data:
# Fixed holdback samples for validation:
holdback.samples <- sample(1:ncol(se),1500)
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

#### Load weights of pretrained model:
cvae %>% load_model_weights_hdf5("TRAINED_MODELS/cvae_deltas_ARCHS4_v212_pretr_512_64_mouse.hdf5")



########################################################################################################  
##########################################    FINE-TUNING      ######################################### 
########################################################################################################  
batch_size <- 512

##### Learning rate scheduling for finetuning: 
lr_schedule <- function(epoch, current_lr) {
    if (epoch <= 100) lr <- 1e-5
    else {lr <- 1e-6} #First LR drop
    return(lr)
}

lr_sch <- callback_learning_rate_scheduler(lr_schedule)

##### Early stopping callback:
early_stopping <- callback_early_stopping(monitor = "val_loss", min_delta = 0,
                                          patience = 100, verbose = 0, mode = "auto",
                                          baseline = NULL, restore_best_weights = TRUE)

#### Checkpoint callback
checkpoint <- callback_model_checkpoint(
    "TRAINED_MODELS/cvae_deltas_ARCHS4_v212_FT_512_64_mouse.hdf5",
    monitor = "val_loss",
    verbose = 0,
    save_best_only = TRUE,
    save_weights_only = TRUE,
    mode = c( "min"),
    period = NULL,
    save_freq = "epoch"
)

###### Training:
nepochs <- 200

######  Fit the model using also our specified callbacks for scheduling and early stopping:
history <- cvae %>% fit(
    x=list(delta_input=train_xD,CONTEXT=train_xC),
    y=train_xD, 
    shuffle = TRUE, 
    epochs = nepochs,
    batch_size = batch_size, 
    validation_data=list(list(test_xD,test_xC),test_xD), 
    callbacks = list( early_stopping, lr_sch, checkpoint )
)

#saveRDS(history,"TRAINED_MODELS/cvae_deltas_ARCHS4_v212_FT_512_64_mouse.rds")
cvae %>% load_model_weights_hdf5("TRAINED_MODELS/cvae_deltas_ARCHS4_v212_FT_512_64_mouse.hdf5")
save_model_hdf5(encoderD,"TRAINED_MODELS/DeltaEncoder_ARCHS_v212_64_mouse.hdf5")
save_model_hdf5(generatorD,"TRAINED_MODELS/DeltaDecoder_ARCHS_v212_64_mouse.hdf5")





















