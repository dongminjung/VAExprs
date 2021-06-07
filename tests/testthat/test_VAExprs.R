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
mat[mat < 0] <- 0
x_train <- as.matrix(t(mat)) 

batch_size <- 32
original_dim <- 1000
intermediate_dim <- 512
epochs <- 2
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
                      use_generator = FALSE)

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
                                    use_generator = FALSE)

gen_sample_result <- gen_exprs(vae_result, num_samples = 100)



test_that("fit_vae: fit_vae yields a model", {
  expect_type(vae_result$model, "closure")
})



test_that("fit_vae: fit_vae yields an encoder", {
  expect_type(vae_result$encoder, "closure")
})



test_that("fit_vae: fit_vae yields a decoder", {
  expect_type(vae_result$decoder, "closure")
})



test_that("fit_vae: fit_vae from preprocessed result", {
  expect_type(vae_result_preprocessing, "list")
})



test_that("fit_vae: add layer for encoder", {
  expect_type(fit_vae(x_train = x_train,
                      encoder_layers = list(layer_input(shape = c(original_dim)),
                                            layer_dense(units = intermediate_dim,
                                                        activation = "relu"),
                                            layer_dense(units = intermediate_dim,
                                                        activation = "relu")),
                      decoder_layers = list(layer_dense(units = intermediate_dim,
                                                        activation = "relu"),
                                            layer_dense(units = original_dim,
                                                        activation = "sigmoid")),
                      epochs = epochs, batch_size = batch_size,
                      validation_split = 0.5,
                      use_generator = FALSE), "list")
})



test_that("fit_vae: dimension of latent vector", {
  expect_type(fit_vae(x_train = x_train,
                      encoder_layers = list(layer_input(shape = c(original_dim)),
                                            layer_dense(units = intermediate_dim,
                                                        activation = "relu"),
                                            layer_dense(units = intermediate_dim,
                                                        activation = "relu")),
                      decoder_layers = list(layer_dense(units = intermediate_dim,
                                                        activation = "relu"),
                                            layer_dense(units = original_dim,
                                                        activation = "sigmoid")),
                      latent_dim = 3,
                      epochs = epochs, batch_size = batch_size,
                      validation_split = 0.5,
                      use_generator = FALSE), "list")
})



test_that("fit_vae: miss an encoder", {
  expect_error(
    fit_vae(x_train = x_train,
            decoder_layers = list(layer_dense(units = intermediate_dim,
                                              activation = "relu"),
                                  layer_dense(units = original_dim,
                                              activation = "sigmoid")),
            epochs = epochs, batch_size = batch_size,
            validation_split = 0.5,
            use_generator = FALSE)
  )
})



test_that("fit_vae: negative regularization parameter", {
  expect_error(
    fit_vae(x_train = x_train,
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
            regularization = -1)
  )
})



test_that("fit_vae: miss epochs", {
  expect_error(
    fit_vae(x_train = x_train,
            encoder_layers = list(layer_input(shape = c(original_dim)),
                                  layer_dense(units = intermediate_dim,
                                              activation = "relu")),
            decoder_layers = list(layer_dense(units = intermediate_dim,
                                              activation = "relu"),
                                  layer_dense(units = original_dim,
                                              activation = "sigmoid")),
            batch_size = batch_size,
            validation_split = 0.5,
            use_generator = FALSE)
  )
})



test_that("fit_vae: miss batch_size", {
  expect_error(
    fit_vae(x_train = x_train,
            encoder_layers = list(layer_input(shape = c(original_dim)),
                                  layer_dense(units = intermediate_dim,
                                              activation = "relu")),
            decoder_layers = list(layer_dense(units = intermediate_dim,
                                              activation = "relu"),
                                  layer_dense(units = original_dim,
                                              activation = "sigmoid")),
            epochs = epochs,
            validation_split = 0.5,
            use_generator = FALSE)
  )
})



test_that("plot_vae: output", {
  expect_visible(plot_vae(vae_result$model))
})



test_that("plot_vae: model is required", {
  expect_error(plot_model(vae_result))
})



test_that("gen_exprs: result", {
  expect_type(gen_sample_result, "list")
})



test_that("gen_exprs: miss number of samples", {
  expect_error(gen_exprs(vae_result))
})



test_that("plot_aug: output", {
  expect_visible(plot_aug(gen_sample_result, "PCA"))
})



test_that("plot_aug: miss plot_fun", {
  expect_error(plot_aug(gen_sample_result))
})
