\name{gen_exprs}
\alias{gen_exprs}

\title{
Generate samples with expression data
}

\description{
This function generate expression data by drawing samples from the latent vectors following the standard multivariate Gaussian distribution (the standard multivariate normal distribution) for convenience. However, this assumption for the prior may not be appropriate because there may be underlying distinctions between groups of samples. Any density function can be modeled by the Gaussian mixture model. Here, by using the library "mclust", the finite Gaussian mixture is applied for such sampling. Note that the Gaussian mixture model is not used for fitting in the function "fit_vae".
}

\usage{
gen_exprs(x, num_samples,
        batch_size, use_generator = FALSE)
}

\arguments{
    \item{x}{
result of the function "fit_vae"
}
    \item{num_samples}{
number of samples to be generated
}
    \item{batch_size}{
batch size
}
    \item{use_generator}{
use data generator if TRUE (default: FALSE)
}
}


\value{
    \item{x_gen}{
generated expression data, where each row is a cell and each column is a gene
}
    \item{y_gen}{
geneated labels
}
    \item{x_train}{
real expression data, where each row is a cell and each column is a gene
}
    \item{y_train}{
real labels
}
    \item{latent_vector}{
latent vector from real expression data
}
}


\author{
Dongmin Jung
}


\seealso{
mclust::mclustBIC, mclust::mclustModel, mclust::sim, DeepPINCS::multiple_sampling_generator, gradDescent::minmaxDescaling, CatEncoders::inverse.transform
}


\examples{
### simulate differentially expressed genes
set.seed(1)
g <- 3
n <- 100
m <- 1000
mu <- 5
sigma <- 5
mat <- matrix(rnorm(n*m*g, mu, sigma), m, n*g)
rownames(mat) <- paste0("gene", seq_len(m))
colnames(mat) <- paste0("cell", seq_len(n*g))
group <- factor(sapply(seq_len(g), function(x) { 
    rep(paste0("group", x), n)
}))
names(group) <- colnames(mat)
mu_upreg <- 6
sigma_upreg <- 10
deg <- 100
for (i in seq_len(g)) {
    mat[(deg*(i-1) + 1):(deg*i), group == paste0("group", i)] <- 
        mat[1:deg, group==paste0("group", i)] + rnorm(deg, mu_upreg, sigma_upreg)
}
# positive expression only
mat[mat < 0] <- 0
x_train <- as.matrix(t(mat))


### model
batch_size <- 32
original_dim <- 1000
intermediate_dim <- 512
epochs <- 2
# VAE
vae_result <- fit_vae(x_train = x_train,
                    encoder_layers = list(layer_input(shape = c(original_dim)),
                                        layer_dense(units = intermediate_dim,
                                                    activation = "relu")),
                    decoder_layers = list(layer_dense(units = intermediate_dim,
                                                    activation = "relu"),
                                        layer_dense(units = original_dim,
                                                    activation = "sigmoid")),
                    epochs = epochs, batch_size = batch_size,
                    validation_split = 0.5,
                    use_generator = FALSE,
                    callbacks = keras::callback_early_stopping(
                        monitor = "val_loss",
                        patience = 10,
                        restore_best_weights = TRUE))
# plot
plot_vae(vae_result$model)


### generate samples
set.seed(1)
gen_sample_result <- gen_exprs(vae_result, num_samples = 100)
}
