reticulate::use_virtualenv("/tungstenfs/groups/gbioinfo/sharedSoft/virtualenvs/r-4.1-bioc-3.13-reticulate-keras-2.4.0-tensorflow-2.4.0-gpu/")
Sys.setenv("CUDA_VISIBLE_DEVICES" = "0" ) # Define visible  GPU devices
ngpus=length(strsplit(Sys.getenv("CUDA_VISIBLE_DEVICES"),",")[[1]])
reticulate::py_config()
library(keras)
K <- backend() # manual add-on
library(tensorflow)
library(Matrix)
library(SummarizedExperiment)


##########################################
###### Data loading and preprocessing ####
##########################################

### Read in summarized experiment:
se <-readRDS("mouse_matrix_v212_237Kx20K_se.rds")

### Create arrays for context (lcpm) data :
temp <- assays(se)[[1]]
pc <- 4
x <- sweep(temp, 2, colSums(temp), FUN = "/") * 1e+06
x <- log2(x+pc)
x <- t(x)                                 #Need to transpose before passing to the model
rm(temp)

####### Splitting in training and validation data:
holdback.samples <- sample(1:nrow(se),8000) 
train.se <- se[,-holdback.samples]
test.se  <- se[, holdback.samples]
rm(se)
train_x <- x[-holdback.samples,]
test_x  <- x[ holdback.samples,]
rm(x)
gc()



##########################################
#### Model definition and compilation ####
##########################################

### Load context VAE model:
ngenes <- ncol(train_x)
if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
source(context_VAE.R)


### VAE loss function:
vae_loss <- function(x, x_decoded_mean){
    reconstruction_loss  <-  loss_mean_squared_error(x, x_decoded_mean)
    kl_loss <- -kl_weight*0.5*K$mean(1 + z_log_var-log_var_prior - K$square(z_mean)/var_prior - K$exp(z_log_var)/var_prior, axis = -1L)  # More general formula
    reconstruction_loss + kl_loss
}


#### Custom correlation function to keep track of:
cor_metric <- function(y_true, y_pred) {  
    x = y_true
    y = y_pred
    xm = x-K$mean(x)
    ym = y-K$mean(y)
    r_num = K$sum(tf$multiply(xm,ym))
    r_den = K$sqrt(tf$multiply(K$sum(K$square(xm)), K$sum(K$square(ym))))
    r = r_num / r_den
    r = K$maximum(K$minimum(r, 1.0), -1.0)
    return (K$square(r))
}


#### Compile model
vae %>% compile(
    loss = vae_loss,
    optimizer = "adam",
    metrics = custom_metric("cor",cor_metric)
)




########################################################################################################  
#######################################     TRAINING REGIME      ####################################### 
########################################################################################################  

##### Learning rate scheduler: 
burn_in.nepochs <- 40 
burn_in_lr <- 1e-5  
batch_size <- 512

#### LR scheduler:
lr_schedule <- function(epoch, current_lr) {
    if (epoch <= burn_in.nepochs ) lr <- burn_in_lr
    else if (epoch < 250) lr <- min (5e-4, burn_in_lr + 1.5e-4 * ( (epoch-burn_in.nepochs)/20) )#Increase lr linearly up to 5e-4
    else if (epoch > 300 ) {lr <- 3e-6 } #Second LR drop  (Cool down)
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
    "TRAINED_MODELS/vae_lcpm_ARCHS_v212_64_mouse.hdf5",
    monitor = "val_loss",
    verbose = 0,
    save_best_only = TRUE,
    save_weights_only = TRUE,
    mode = c( "min"),
    period = NULL,
    save_freq = "epoch"
)

###### Training:
nepochs <- 400
history <- vae %>% fit(
    x=train_x,
    y=train_x, 
    shuffle = TRUE, 
    epochs = nepochs,
    batch_size = batch_size, 
    validation_data=list(test_x,test_x), 
    callbacks = list( early_stopping, lr_sch,checkpoint )
)

#saveRDS(history,"Trained_models/history_vae_deJUNKER_lcpm_ARCHS_v212_64_mouse.rds")
vae %>% load_model_weights_hdf5("TRAINED_MODELS/vae_lcpm_ARCHS_v212_64_mouse.hdf5")
save_model_hdf5(encoder,"TRAINED_MODELS/ContextEncoder_ARCHS_v212_64_mouse.hdf5")








