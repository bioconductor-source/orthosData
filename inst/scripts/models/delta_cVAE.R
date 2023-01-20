##### conditionalVAE model for  contrast data (delta lcpms)
##### Layer width: 20K+64 | 2048 | 1024 |--512+64--| 1024 | 2048 | 20K

# Model parameters 
latentD <- 512L # Latent dimension for Delta encoder
latentC <- 64L # Latent dimension for Context encoder
drop_rate <- 0.2 # Dropout rate
gene_dim <- ngenes # Number of features (genes) in the dataset
epsilon_std <- 0.5 # Standard deviation of the prior latent distribution (vanilla = 1)
var_prior <- epsilon_std**2
log_var_prior <- log(var_prior)
kl_weight <- 0.20 # Weight for the Delta Kulllback-Leibler divergence loss (vanilla = 1)

### Define encoder input layers:
xD <- keras::layer_input(shape = c(gene_dim), name = "delta_input")
cont_input <- keras::layer_input(shape = c(latentC), name = 'CONTEXT')
all_inputs <- keras::layer_concatenate(list(xD, cont_input))

#### Delta encoder definition:
hD <- keras::layer_dense(all_inputs, 4 * latentD, activation = "elu")
hD <- keras::layer_dropout(hD, rate = drop_rate)
hD <- keras::layer_dense(hD, 2 * latentD, activation = "elu")
hD <- keras::layer_dropout(hD, rate = drop_rate)
z_meanD <- keras::layer_dense(hD, latentD)
z_log_varD <- keras::layer_dense(hD, latentD)

# Define delta encoder
encoderD <- keras::keras_model(inputs = c(xD, cont_input), outputs = z_meanD)

#### Sampling from the Delta latent space:
samplingD <- function(arg) {
    z_meanD <- arg[, seq_len(latentD)]
    z_log_varD <- arg[, (latentD + 1):(2 * latentD)]
    epsilonD <- K$random_normal(
        shape = c(K$shape(z_meanD)[[1]]),
        mean = 0.,
        stddev = epsilon_std
    )
    z_meanD + K$exp(z_log_varD/2) * epsilonD
}

# Lambda layer for variational sampling:
zD <- keras::layer_concatenate(list(z_meanD, z_log_varD)) %>%
    keras::layer_lambda(samplingD)

# Merge delta latent space with context latent space:
z_concat <- keras::layer_concatenate(list(zD, cont_input))

# Define layers for the Delta decoder (no batch-norm. Seems to increase overfitting):
decoder_h1 <- keras::layer_dense(units = 2*latentD, activation = "elu")
decoder_h2 <- keras::layer_dropout(rate = drop_rate)
decoder_h3 <- keras::layer_dense(units = 4*latentD, activation = "elu")
decoder_h4 <- keras::layer_dropout(rate = drop_rate)

decoder_out <- keras::layer_dense(units = ngenes, activation = "linear")

h_p <- decoder_h1(z_concat)
h_p <- decoder_h2(h_p)
h_p <- decoder_h3(h_p)
h_p <- decoder_h4(h_p)
outputs <- decoder_out(h_p)

# Define full cvae (Input: lcpm, deltas  Target Outputs: Decoded deltas)
cvae <- keras::keras_model(inputs = c(xD, cont_input), outputs)

# Reuse decoder layers to define the generator (concatenated latent space to gene output) separately
d_in <- layer_input(shape = latentC + latentD)
d_h <- decoder_h1(d_in)
d_h <- decoder_h2(d_h)
d_h <- decoder_h3(d_h)
d_h <- decoder_h4(d_h)
d_out <- decoder_out(d_h)
generatorD <- keras::keras_model(d_in, d_out)
