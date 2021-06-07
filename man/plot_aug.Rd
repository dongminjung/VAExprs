\name{plot_aug}
\alias{plot_aug}

\title{
Visualization for augmented data
}

\description{
For augmented data, we can create plots for specific types of dimension reduction.
}

\usage{
plot_aug(x, plot_fun, ...)
}

\arguments{
    \item{x}{
result of the function "gen_exprs"
}
    \item{plot_fun}{
"PCA", "MDS", "TSNE", "UMAP", "NMF", or "DiffusionMap"
}
    \item{...}{
additional parameters for the reduced dimension plots such as "scater::runPCA"
}
}


\value{
plot for augmented data
}


\author{
Dongmin Jung
}


\seealso{
SingleCellExperiment::SingleCellExperiment, scater::logNormCounts, scater::runPCA, scater::runMDS, scater::runTSNE, scater::runUMAP, scater::runNMF, scater::runDiffusionMap, scater::plotPCA, scater::plotMDS, scater::plotTSNE, scater::plotUMAP, scater::plotNMF, scater::plotDiffusionMap
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
# plot
plot_aug(gen_sample_result, "PCA")
}
