# VAExprs

## Introduction

A fundamental problem in biomedical research is the low number of observations, mostly due to a lack of available biosamples, prohibitive costs, or ethical reasons. By augmenting a few real observations with artificially generated samples, their analysis could lead to more robust and higher reproducible. One possible solution to the problem is the use of generative models, which are statistical models of data that attempt to capture the entire probability distribution from the observations. Using the variational autoencoder (VAE), a well-known deep generative model, this package is aimed to generate samples with gene expression data, especially for single-cell RNA-seq data. Furthermore, the VAE can use conditioning to produce specific cell types or subpopulations. The conditional VAE (CVAE) allows us to create targeted samples rather than completely random ones.





## Installation

Before installing this package, tensorflow and kears must be installed in Python and connected to R (https://tensorflow.rstudio.com, https://keras.rstudio.com).
To install the development version from GitHub:
``` 
devtools::install_github("dongminjung/VAExprs")
```





## Tutorial

The "yan" data set is single-cell RNA sequencing data with 20214 genes and 90 cells from human preimplantation embryos and embryonic stem cells at different passages. The rows in the dataset correspond to genes and columns correspond to cells. The "SingleCellExperiment" class can be used to store and manipulate single-cell genomics data. It extends the "RangedSummarizedExperiment" class and follows similar conventions. The object "sce" can be created by the data "yan" with cell type annotation "ann".

```
library(VAExprs)
library(SC3)
library(SingleCellExperiment)

# create a SingleCellExperiment object
sce <- SingleCellExperiment(
    assays = list(counts = as.matrix(yan)),
    colData = ann
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
# remove genes that are not expressed in any samples
sce <- sce[which(rowMeans(assay(sce)) > 0),]
dim(assay(sce))
```

The CVAE model can be built by using the function "fit_vae" with gene expression data and the cell annotation from the object "sce". The overall loss function of the VAE is expressed as a weighted sum of the reconstruction loss and the regularization loss. The reconstruction loss is the binary cross-entropy loss between the input and output and the regularization loss is simply the Kullback-Leibler divergence measure.
```
# model parameters
batch_size <- 32
original_dim <- 19595
intermediate_dim <- 256
epochs <- 100

# model
cvae_result <- fit_vae(object = sce,
                       encoder_layers = list(layer_input(shape = c(original_dim)),
                                             layer_dense(units = intermediate_dim,
                                                         activation = "relu")),
                       decoder_layers = list(layer_dense(units = intermediate_dim,
                                                         activation = "relu"),
                                             layer_dense(units = original_dim,
                                                         activation = "sigmoid")),
                       epochs = epochs, batch_size = batch_size,
                       use_generator = TRUE,
                       callbacks = keras::callback_early_stopping(
                           monitor = "loss",
                           patience = 20,
                           restore_best_weights = TRUE))

Epoch 1/100
3/3 [==============================] - 2s 523ms/step - loss: 13073.0146
Epoch 2/100
3/3 [==============================] - 1s 221ms/step - loss: 10780.1800
Epoch 3/100
3/3 [==============================] - 1s 213ms/step - loss: 9844.0072
Epoch 4/100
3/3 [==============================] - 1s 243ms/step - loss: 9296.4753
Epoch 5/100
3/3 [==============================] - 1s 198ms/step - loss: 8945.9284
Epoch 6/100
3/3 [==============================] - 1s 235ms/step - loss: 8680.7949
Epoch 7/100
3/3 [==============================] - 1s 230ms/step - loss: 8421.7432
Epoch 8/100
3/3 [==============================] - 1s 201ms/step - loss: 8242.2995
Epoch 9/100
3/3 [==============================] - 1s 194ms/step - loss: 8718.2835
Epoch 10/100
3/3 [==============================] - 1s 229ms/step - loss: 8247.0957
Epoch 11/100
3/3 [==============================] - 1s 199ms/step - loss: 8338.7660
Epoch 12/100
3/3 [==============================] - 1s 206ms/step - loss: 8325.5381
Epoch 13/100
3/3 [==============================] - 1s 228ms/step - loss: 8055.5200
Epoch 14/100
3/3 [==============================] - 1s 209ms/step - loss: 8247.4580
Epoch 15/100
3/3 [==============================] - 1s 392ms/step - loss: 8028.6133
Epoch 16/100
3/3 [==============================] - 1s 218ms/step - loss: 8177.4984
Epoch 17/100
3/3 [==============================] - 1s 205ms/step - loss: 7902.1089
Epoch 18/100
3/3 [==============================] - 1s 194ms/step - loss: 8087.9927
Epoch 19/100
3/3 [==============================] - 1s 192ms/step - loss: 7990.8016
Epoch 20/100
3/3 [==============================] - 1s 219ms/step - loss: 7896.8932
Epoch 21/100
3/3 [==============================] - 1s 193ms/step - loss: 7941.2196
Epoch 22/100
3/3 [==============================] - 1s 201ms/step - loss: 7940.5312
Epoch 23/100
3/3 [==============================] - 1s 202ms/step - loss: 7982.6789
Epoch 24/100
3/3 [==============================] - 1s 201ms/step - loss: 7966.7583
Epoch 25/100
3/3 [==============================] - 1s 209ms/step - loss: 8028.1893
Epoch 26/100
3/3 [==============================] - 1s 208ms/step - loss: 7826.7415
Epoch 27/100
3/3 [==============================] - 1s 375ms/step - loss: 7874.1709
Epoch 28/100
3/3 [==============================] - 1s 202ms/step - loss: 7739.7087
Epoch 29/100
3/3 [==============================] - 1s 205ms/step - loss: 7717.1283
Epoch 30/100
3/3 [==============================] - 1s 216ms/step - loss: 7767.9258
Epoch 31/100
3/3 [==============================] - 1s 198ms/step - loss: 7774.4969
Epoch 32/100
3/3 [==============================] - 1s 203ms/step - loss: 7796.2261
Epoch 33/100
3/3 [==============================] - 1s 214ms/step - loss: 7869.4207
Epoch 34/100
3/3 [==============================] - 1s 379ms/step - loss: 7911.8880
Epoch 35/100
3/3 [==============================] - 1s 198ms/step - loss: 7946.8022
Epoch 36/100
3/3 [==============================] - 1s 195ms/step - loss: 7778.8597
Epoch 37/100
3/3 [==============================] - 1s 199ms/step - loss: 7760.7996
Epoch 38/100
3/3 [==============================] - 1s 217ms/step - loss: 7835.4167
Epoch 39/100
3/3 [==============================] - 1s 210ms/step - loss: 7901.9056
Epoch 40/100
3/3 [==============================] - 1s 210ms/step - loss: 7803.0190
Epoch 41/100
3/3 [==============================] - 1s 197ms/step - loss: 7802.9461
Epoch 42/100
3/3 [==============================] - 1s 197ms/step - loss: 7717.3841
Epoch 43/100
3/3 [==============================] - 1s 196ms/step - loss: 7753.0926
Epoch 44/100
3/3 [==============================] - 1s 197ms/step - loss: 7818.5275
Epoch 45/100
3/3 [==============================] - 1s 198ms/step - loss: 7731.0760
Epoch 46/100
3/3 [==============================] - 1s 224ms/step - loss: 7658.5674
Epoch 47/100
3/3 [==============================] - 1s 220ms/step - loss: 7693.5033
Epoch 48/100
3/3 [==============================] - 1s 200ms/step - loss: 7844.3693
Epoch 49/100
3/3 [==============================] - 1s 200ms/step - loss: 7780.8586
Epoch 50/100
3/3 [==============================] - 1s 201ms/step - loss: 7702.2021
Epoch 51/100
3/3 [==============================] - 1s 206ms/step - loss: 7785.0112
Epoch 52/100
3/3 [==============================] - 1s 219ms/step - loss: 7800.9276
Epoch 53/100
3/3 [==============================] - 1s 200ms/step - loss: 7751.7611
Epoch 54/100
3/3 [==============================] - 1s 195ms/step - loss: 7716.8810
Epoch 55/100
3/3 [==============================] - 1s 196ms/step - loss: 7704.1950
Epoch 56/100
3/3 [==============================] - 1s 200ms/step - loss: 7755.0158
Epoch 57/100
3/3 [==============================] - 1s 215ms/step - loss: 7528.5760
Epoch 58/100
3/3 [==============================] - 1s 214ms/step - loss: 7717.8184
Epoch 59/100
3/3 [==============================] - 1s 201ms/step - loss: 7711.2790
Epoch 60/100
3/3 [==============================] - 1s 193ms/step - loss: 7721.7414
Epoch 61/100
3/3 [==============================] - 1s 199ms/step - loss: 7849.1649
Epoch 62/100
3/3 [==============================] - 1s 202ms/step - loss: 7674.6919
Epoch 63/100
3/3 [==============================] - 1s 386ms/step - loss: 7544.1912
Epoch 64/100
3/3 [==============================] - 1s 201ms/step - loss: 7586.4850
Epoch 65/100
3/3 [==============================] - 1s 201ms/step - loss: 7546.8431
Epoch 66/100
3/3 [==============================] - 1s 201ms/step - loss: 7567.8174
Epoch 67/100
3/3 [==============================] - 1s 232ms/step - loss: 7476.5628
Epoch 68/100
3/3 [==============================] - 1s 211ms/step - loss: 7665.1159
Epoch 69/100
3/3 [==============================] - 1s 209ms/step - loss: 7533.3299
Epoch 70/100
3/3 [==============================] - 1s 212ms/step - loss: 7611.5394
Epoch 71/100
3/3 [==============================] - 1s 198ms/step - loss: 7685.3724
Epoch 72/100
3/3 [==============================] - 1s 199ms/step - loss: 7634.1060
Epoch 73/100
3/3 [==============================] - 1s 207ms/step - loss: 7655.5771
Epoch 74/100
3/3 [==============================] - 1s 201ms/step - loss: 7682.9412
Epoch 75/100
3/3 [==============================] - 1s 194ms/step - loss: 7784.4570
Epoch 76/100
3/3 [==============================] - 1s 223ms/step - loss: 7652.0124
Epoch 77/100
3/3 [==============================] - 1s 195ms/step - loss: 7527.1865
Epoch 78/100
3/3 [==============================] - 1s 195ms/step - loss: 7698.5597
Epoch 79/100
3/3 [==============================] - 1s 199ms/step - loss: 7602.8826
Epoch 80/100
3/3 [==============================] - 1s 203ms/step - loss: 7746.2612
Epoch 81/100
3/3 [==============================] - 1s 208ms/step - loss: 7578.2534
Epoch 82/100
3/3 [==============================] - 1s 206ms/step - loss: 7525.1266
Epoch 83/100
3/3 [==============================] - 1s 199ms/step - loss: 7421.0713
Epoch 84/100
3/3 [==============================] - 1s 195ms/step - loss: 7575.2210
Epoch 85/100
3/3 [==============================] - 1s 200ms/step - loss: 7645.3781
Epoch 86/100
3/3 [==============================] - 1s 211ms/step - loss: 7578.3859
Epoch 87/100
3/3 [==============================] - 1s 222ms/step - loss: 7730.1956
Epoch 88/100
3/3 [==============================] - 1s 199ms/step - loss: 7755.5093
Epoch 89/100
3/3 [==============================] - 1s 198ms/step - loss: 7631.6792
Epoch 90/100
3/3 [==============================] - 1s 190ms/step - loss: 7554.1834
Epoch 91/100
3/3 [==============================] - 1s 211ms/step - loss: 7639.7676
Epoch 92/100
3/3 [==============================] - 1s 221ms/step - loss: 7534.7118
Epoch 93/100
3/3 [==============================] - 1s 195ms/step - loss: 7695.3576
Epoch 94/100
3/3 [==============================] - 1s 198ms/step - loss: 7604.4395
Epoch 95/100
3/3 [==============================] - 1s 203ms/step - loss: 7416.3944
Epoch 96/100
3/3 [==============================] - 1s 207ms/step - loss: 7577.1263
Epoch 97/100
3/3 [==============================] - 1s 373ms/step - loss: 7643.6406
Epoch 98/100
3/3 [==============================] - 1s 206ms/step - loss: 7426.4823
Epoch 99/100
3/3 [==============================] - 1s 228ms/step - loss: 7471.0667
Epoch 100/100
3/3 [==============================] - 1s 200ms/step - loss: 7493.5422
```

![history plot](docs/history_plot.png)

The function "plot_vae" draws the plot for model architecture.
```
# model architecture
plot_vae(cvae_result$model)
```

![history plot](docs/model_architecture.png)

The function "gen_exprs" can generate samples with expression data by using the trained model.
```
# sample generation
set.seed(1)
gen_sample_result <- gen_exprs(cvae_result, 100,
                               batch_size, use_generator = TRUE)
```

The function "plot_aug" uses reduced dimension plots for augmented data visualization.
```
# plot for augmented data
plot_aug(gen_sample_result, "PCA")
```

![history plot](docs/augmentation_plot.png)