\name{fit_vae}
\alias{fit_vae}

\title{
Variational autoencoder model fitting
}

\description{
A fundamental problem in biomedical research is the low number of observations available. Augmenting a few real observations with generated in silico samples could lead to more robust analysis. Here, the variational autoencoder (VAE) is used for the realistic generation of single-cell RNA-seq data. Also, the conditional variational autoencoder (CVAE) can be used if labels of samples are available. This function allows us to fit variational autoencoders with the standard Gaussian prior to expression data. It is assumed that there will likely be no clusters in the latent space representation of variational autoencoders.
}

\usage{
fit_vae(object = NULL,
        x_train = NULL,
        x_val = NULL,
        y_train = NULL,
        y_val = NULL,
        encoder_layers,
        decoder_layers,
        latent_dim = 2,
        regularization = 1,
        epochs,
        batch_size,
        preprocessing = list(
            x_train = NULL,
            x_val = NULL,
            y_train = NULL,
            y_val = NULL,
            minmax = NULL,
            lenc = NULL),
        use_generator = FALSE,
        optimizer = "adam",
        validation_split = 0, ...)
}

\arguments{
    \item{object}{
SummarizedExperiment object
}
    \item{x_train}{
expression data for train, where each row is a cell and each column is a gene
}
    \item{x_val}{
expression data for validation, where each row is a cell and each column is a gene
}
    \item{y_train}{
labels for train
}
    \item{y_val}{
labels for validation
}
    \item{encoder_layers}{
list of layers for encoder
}
    \item{decoder_layers}{
list of layers for decoder
}
    \item{latent_dim}{
dimension of latent vector (default: 2)
}
    \item{regularization}{
regularization parameter, which is nonnegative (default: 1)
}
    \item{epochs}{
number of epochs
}
    \item{batch_size}{
batch size
}
    \item{preprocessing}{
list of preprocessed results, they are set to NULL as default
\itemize{
    \item x_train : expression data for train
    \item x_val : expression data for validation
    \item y_train : labels for train
    \item y_val : labels for validation
    \item minmax : result of min-max normalization
    \item lenc : encoded labels
}
}
    \item{use_generator}{
use data generator if TRUE (default: FALSE)
}
    \item{optimizer}{
name of optimizer (default: adam)
}
    \item{validation_split}{
proportion of validation data, it is ignored when there is a validation set (default: 0)
}
    \item{...}{
additional parameters for the "fit" or "fit_generator"
}
}


\value{
    \item{model}{
trained VAE model
}
    \item{encoder}{
trained encoder model
}
    \item{decoder}{
trained decoder model
}
    \item{preprocessing}{
preprocessed results
}
}


\author{
Dongmin Jung
}


\seealso{
SummarizedExperiment::assay, SummarizedExperiment::colData, scater::logNormCounts, gradDescent::minmaxScaling, keras::fit, keras::fit_generator, keras::compile, CatEncoders::LabelEncoder.fit, CatEncoders::transform, DeepPINCS::multiple_sampling_generator
}


\references{
Marouf, M., Machart, P., Bansal, V., Kilian, C., Magruder, D. S., Krebs, C. F., & Bonn, S. (2020). Realistic in silico generation and augmentation of single-cell RNA-seq data using generative adversarial networks. Nature communications, 11(1), 1-12.
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


### from preprocessing
vae_result_preprocessing <- fit_vae(preprocessing = vae_result$preprocessing,
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
}
