\name{plot_vae}
\alias{plot_vae}

\title{
Visualization for the variational autoencoder
}

\description{
You can create a plot of the VAE model. This plot can help you check that the model is connected the way you intended. The node colors indicate the components of the VAE.
}

\usage{
plot_vae(x, node_color = list(encoder_col = "tomato",
                            mean_vector_col = "orange",
                            stddev_vector_col = "lavender",
                            latent_vector_col = "lightblue",
                            decoder_col = "palegreen",
                            condition_col = "gray"))
}

\arguments{
    \item{x}{
VAE model
}
    \item{node_color}{
node colors for encoder(default: tomato), mean vector(default: orange), standard deviation vector(default: lavender), latent_vector(default: lightblue), decoder(default: palegreen), and condition(default: gray)
}
}


\value{
plot for the model architecture
}


\author{
Dongmin Jung
}


\seealso{
purrr::map, purrr::map_chr, purrr::pluck, purrr::imap_dfr, DiagrammeR::grViz
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
}
