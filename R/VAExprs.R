fit_vae <- function(object = NULL,
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
                    validation_split = 0, ...) {
    if (tensorflow::tf$executing_eagerly())
        tensorflow::tf$compat$v1$disable_eager_execution()
    
    result <- NULL
    result$preprocessing <- NULL
    x_val <- NULL
    y_val <- NULL
    if (regularization < 0) {
        stop("regularization parameter should be nonnegative")
    }
    
    
    ### pre-processing
    if (all(unlist(lapply(preprocessing, Negate(is.null))))) {
        result$preprocessing <- preprocessing
        x_train <- preprocessing$x_train
        if (!is.null(preprocessing$x_val)) x_val <- preprocessing$x_val
        if (!is.null(preprocessing$y_train)) y_train <- preprocessing$y_train
        if (!is.null(preprocessing$y_val)) y_val <- preprocessing$y_val
    } else if (!is.null(object)) {
        message("pre-processing...")
        x <- t(SummarizedExperiment::assay(scater::logNormCounts(object), "logcounts"))
        y <- NULL
        if (length(colData(object)) > 0) {
            y <- SummarizedExperiment::colData(object)[,1]
        }
        
        if (validation_split > 0) {
            idx <- sample(seq_len(nrow(x)))
            train_idx <- seq_len(nrow(x)) %in%
                idx[seq_len(round(nrow(x) * (1 - validation_split)))]
            x_train <- x[train_idx,]
            x_val <- x[!train_idx,]
            
            if (!is.null(y)) {
                y_train <- y[train_idx]
                y_val <- y[!train_idx]
            }
        } else {
            x_train <- x
            if (!is.null(y)) {
                y_train <- y
            }
        }
    }
    
    if (any(unlist(lapply(preprocessing, is.null)))) {
        if (is.null(object) & is.null(x_val) & validation_split) {
            message("preprocessing...")
            x <- x_train
            idx <- sample(seq_len(nrow(x)))
            train_idx <- seq_len(nrow(x)) %in%
                idx[seq_len(round(nrow(x) * (1 - validation_split)))]
            x_train <- x[train_idx,]
            x_val <- x[!train_idx,]
            if (!is.null(y_train) & is.null(y_val)) {
                y <- y_train
                y_train <- y[train_idx]
                y_val <- y[!train_idx]
            }
        }
    }
    
    
    ### normalization
    if (any(unlist(lapply(preprocessing, is.null)))) {
        message("normalizing...")
        minmax <- gradDescent::minmaxScaling(data.frame(x_train))
        x_train <- as.matrix(minmax$scaledDataSet)
        if (!is.null(x_val)) {
            x_val <- as.matrix(gradDescent::minmaxScaling(data.frame(x_val))$scaledDataSet)
        }
        result$preprocessing$minmax <- minmax
    }
    
    
    ### building model
    # VAE
    encoded <- NULL
    for (i in seq_len(length(encoder_layers))) {
        if (is.null(encoded)) {
            x <- encoder_layers[[i]]
            if (is.null(y_train)) {
                encoded <- x
            } else {
                condition <- layer_input(shape = c(1), name = "condition")
                encoded <- layer_concatenate(c(x, condition))
            }
        } else {
            encoded <- encoder_layers[[i]](encoded)
        }
    }
    
    assign("z_mean", layer_dense(encoded, latent_dim, name = "z_mean"),
        envir = globalenv())
    assign("z_log_stddev", layer_dense(encoded, latent_dim, name = "z_log_stddev"),
        envir = globalenv())
    
    z <- layer_lambda(list(z_mean, z_log_stddev),
            function(arg) {
                z_mean <- arg[[1]]
                z_log_stddev <- arg[[2]]
                
                epsilon <- k_random_normal(shape = c(k_shape(z_mean)[[1]]),
                                            mean = 0.0, stddev = 1)
                
                z_mean + k_exp(z_log_stddev)*epsilon
            }, name = "latent")
    
    if (is.null(y_train)) {
        decoded <- z
    } else {
        decoded <- layer_concatenate(c(z, condition))
    }
    for (i in seq_len(length(decoder_layers))) {
        decoded <- decoder_layers[[i]](decoded)
    }
    
    if (is.null(y_train)) {
        model <- keras_model(inputs = x, outputs = decoded)
    } else {
        model <- keras_model(inputs = c(x, condition), outputs = decoded)
    }
    
    # encoder
    if (is.null(y_train)) {
        encoder <- keras_model(inputs = x, outputs = z_mean)
    } else {
        encoder <- keras_model(inputs = c(x, condition), outputs = z_mean)
    }
    
    # decoder
    if (is.null(y_train)) {
        decoder_input <- layer_input(shape = latent_dim)
    } else {
        decoder_input <- layer_input(shape = latent_dim + 1)
    }
    decoder_output <- decoder_input
    for (i in seq_len(length(decoder_layers))) {
        decoder_output <- decoder_layers[[i]](decoder_output)
    }
    decoder <- keras_model(decoder_input, decoder_output)
    
    # loss
    vae_loss <- function(x_true, x_pred) {
        xent_loss <- (original_dim / 1.0) * loss_binary_crossentropy(x_true, x_pred)
        kl_loss <- -0.5 * k_mean(1 + z_log_stddev - k_square(z_mean) - k_exp(z_log_stddev), axis = -1L)
        xent_loss + regularization * kl_loss
    }
    
    model %>% keras::compile(optimizer = optimizer, loss = vae_loss)
    
    
    ### training
    message("training...")
    if (is.null(y_train)) {
        # vanilla vae
        if (!use_generator) {
            # without generator
            validation_data <- NULL
            if (!is.null(x_val)) {
                validation_data <- list(x_val, x_val)
            }
            
            model %>% keras::fit(
                x_train, x_train,
                epochs = epochs, 
                batch_size = batch_size, 
                validation_data = validation_data,
                validation_split = validation_split, ...)
        } else {
            # with generator
            validation_data <- NULL
            validation_steps <- NULL
            if (!is.null(x_val)) {
                validation_data <- DeepPINCS::multiple_sampling_generator(
                    list(x_val), x_val,
                    batch_size = batch_size)
                validation_steps <- ceiling(nrow(x_val)/batch_size)
            }
            
            model %>% keras::fit(
                DeepPINCS::multiple_sampling_generator(
                    list(x_train), x_train,
                    batch_size = batch_size),
                steps_per_epoch = ceiling(nrow(x_train)/batch_size),
                epochs = epochs,
                validation_data = validation_data,
                validation_steps = validation_steps, ...)
        }
    } else {
        # conditional vae
        if (any(unlist(lapply(preprocessing, is.null)))) {
            lenc <- CatEncoders::LabelEncoder.fit(y_train)
            y_train <- CatEncoders::transform(lenc, y_train)
            result$preprocessing$lenc <- lenc
        }
        
        if (!use_generator) {
            # without generator
            validation_data <- NULL
            if (!is.null(x_val)) {
                validation_data <- list(list(x_val, y_val), x_val)
            }
            
            model %>% keras::fit(
                list(x_train, y_train), x_train,
                shuffle = TRUE, 
                epochs = epochs, 
                batch_size = batch_size, 
                validation_data = validation_data,
                validation_split = validation_split, ...)
        } else {
            # with generator
            validation_data <- NULL
            validation_steps <- NULL
            if (!is.null(x_val)) {
                validation_data <- DeepPINCS::multiple_sampling_generator(
                    list(x_val, cbind(y_val)), cbind(x_val),
                    batch_size = batch_size)
                validation_steps <- ceiling(nrow(x_val)/batch_size)
            }
            
            model %>% keras::fit(
                DeepPINCS::multiple_sampling_generator(
                    list(x_train, cbind(y_train)), cbind(x_train),
                    batch_size = batch_size),
                steps_per_epoch = ceiling(nrow(x_train)/batch_size),
                epochs = epochs,
                validation_data = validation_data,
                validation_steps = validation_steps, ...)
        }
    }
    
    rm(z_mean, envir = globalenv())
    rm(z_log_stddev, envir = globalenv())
    
    result$model <- model
    result$encoder <- encoder
    result$decoder <- decoder
    result$preprocessing$x_train <- x_train
    if (!is.null(y_train)) result$preprocessing$y_train <- y_train
    if (!is.null(x_val)) result$preprocessing$x_val <- x_val
    if (!is.null(y_val)) result$preprocessing$y_val <- y_val
    result
}





plot_vae <- function(x, node_color = list(encoder_col = "tomato",
                                        mean_vector_col = "orange",
                                        stddev_vector_col = "lavender",
                                        latent_vector_col = "lightblue",
                                        decoder_col = "palegreen",
                                        condition_col = "gray")) {
    # node information
    layer_name <- x$layers %>%
        purrr::map_chr(~(purrr::`%||%`(purrr::pluck(., "name"), "")))
    layer_type <- purrr::map(x$layers, c(class, 1)) %>%
        unlist() %>%
        strsplit(split = "\\.") %>%
        purrr::map(c(rev, 1)) %>%
        unlist()
    layer_input_shape <- unlist(lapply(x$layers, function(x) {
        ifelse(length(purrr::pluck(x$input_shape, 1)) > 0,
            paste("[", paste(unlist(lapply(x$input_shape,
                function(x) paste("(", toString(paste(x)), ")", sep = ""))),
                    collapse = ", "), "]", sep = ""),
            paste("(", toString(paste(x$input_shape)), ")", sep = ""))
    }))
    layer_input_shape <- gsub("NULL", "None", layer_input_shape)
    layer_output_shape <- unlist(lapply(x$layers, function(x) {
        ifelse(length(purrr::pluck(x$output_shape, 1)) > 0,
            paste("[", paste(unlist(lapply(x$output_shape,
                function(x) paste("(", toString(paste(x)), ")", sep = ""))),
                    collapse = ", "), "]", sep = ""),
            paste("(", toString(paste(x$output_shape)), ")", sep = ""))
    }))
    layer_output_shape <- gsub("NULL", "None", layer_output_shape)
    layer_color <- c(rep(node_color[["encoder_col"]], which(layer_name == "z_mean") - 1),
                    node_color[["mean_vector_col"]],
                    node_color[["stddev_vector_col"]],
                    node_color[["latent_vector_col"]],
                    rep(node_color[["decoder_col"]],
                        length(x$layers) - which(layer_name == "latent")))
    if (any(layer_name == "condition")) {
        layer_color[which(layer_name == "condition")] <- node_color[["condition_col"]]
    }
    node_info <- data.frame(layer_name, layer_type, layer_input_shape, layer_output_shape)
    
    # edge information
    inbound <- lapply(x$layers %>%
                        purrr::map("inbound_nodes") %>%
                        purrr::map(1) %>%
                        purrr::map("inbound_layers"),
                    function(x) switch(length(x), x$name, unlist(x %>% purrr::map("name"))))
    
    if (length(Filter(Negate(is.null), inbound)) == 0) {
        edge_info <- embed(rownames(node_info), dimension = 2)
        from <- edge_info[, 2]
        to <- edge_info[, 1]
        edge_info <- data.frame(from, to, stringsAsFactors = FALSE)
    } else {
        names(inbound) <- layer_name
        inbound <- Filter(Negate(is.null), inbound)
        edge_info <- purrr::imap_dfr(inbound, ~data.frame(from = .x, to = .y, stringsAsFactors = FALSE))
        edge_info$from <- rownames(node_info)[match(edge_info$from, node_info$layer_name)]
        edge_info$to <- rownames(node_info)[match(edge_info$to, node_info$layer_name)]
    }
    
    # plot
    nodes <- paste(unlist(lapply(seq_len(nrow(node_info)),
                                function(x) paste(toString(x),
                                                " [label = '@@",
                                                toString(x),
                                                "', fillcolor = ",
                                                layer_color[x], "]",
                                                sep = ""))), collapse = "")
    names <- paste(unlist(lapply(seq_len(nrow(node_info)), function(x) {
        paste(" [", toString(x), "]: ", "'", node_info$layer_name[x], " ",
            ifelse(node_info$layer_name_sub[x] != "",
                paste("(", node_info$layer_name_sub[x], ")", sep = ""), ""),
            " : ",
            node_info$layer_type[x], " ",
            ifelse(node_info$layer_type_sub[x] != "",
                paste("(", node_info$layer_type_sub[x], ")", sep = ""), ""),
            "|{input: | output:}", "|{", 
            node_info$layer_input_shape[x], "|", node_info$layer_output_shape[x], "}", "'", sep = "")
    })), collapse = "\n")
    edges <- gsub(",", "->", paste(apply(edge_info, 1, toString), collapse = " "))
    DiagrammeR::grViz(paste("digraph{
                            graph [layout = dot]
                            node [shape = record, color = black, style = filled]", 
                            nodes, edges, "} \n", names))
}





gen_exprs <- function(x,
                    num_samples,
                    batch_size,
                    use_generator = FALSE) {
    result <- NULL
    encoder <- x$encoder
    decoder <- x$decoder
    x_train <- x$preprocessing$x_train
    y_train <- x$preprocessing$y_train
    minmax <- x$preprocessing$minmax
    
    x_train_encoded <- predict(encoder, list(x_train, y_train))
    rownames(x_train_encoded) <- rownames(x_train)
    colnames(x_train_encoded) <- paste("latent", seq_len(ncol(x_train_encoded)))
    
    
    ### generating
    message("generating...")
    if (is.null(y_train)) {
        # vanilla vae
        y_gen <- NULL
        model_BIC <- mclust::mclustBIC(x_train_encoded, verbose = FALSE)
        mod <- mclust::mclustModel(x_train_encoded, model_BIC)
        z_sample <- mclust::sim(modelName = mod$modelName, 
                                parameters = mod$parameters, 
                                n = num_samples)
        
        if (!use_generator) {
            x_gen <- predict(decoder, z_sample[,-1])
        } else {
            x_gen <- predict_generator(decoder,
                                        DeepPINCS::multiple_sampling_generator(
                                            list(z_sample[,-1]),
                                            batch_size = batch_size,
                                            shuffle = FALSE),
                                        steps = ceiling(nrow(z_sample)/batch_size))
        }
        
    } else {
        # conditional vae
        all_labels <- names(table(y_train))
        y_gen <- sample(all_labels, num_samples, replace = TRUE,
                        prob = as.vector(table(y_train)/length(y_train)))
        z_sample <- matrix(0, num_samples, ncol(x_train_encoded))
        for (i in seq_len(length(all_labels))) {
            model_BIC <- mclust::mclustBIC(x_train_encoded[y_train == all_labels[i],],
                                        verbose = FALSE)
            mod <- mclust::mclustModel(x_train_encoded[y_train == all_labels[i],], model_BIC)
            temp_z <- mclust::sim(modelName = mod$modelName, 
                                parameters = mod$parameters, 
                                n = table(y_gen)[i])
            z_sample[y_gen == all_labels[i],] <- temp_z[,-1]
        }
        
        if (!use_generator) {
            x_gen <- predict(decoder, cbind(z_sample, y_gen))
        } else {
            x_gen <- predict_generator(decoder,
                                        DeepPINCS::multiple_sampling_generator(
                                            list(cbind(z_sample, y_gen)),
                                            batch_size = batch_size,
                                            shuffle = FALSE),
                                        steps = ceiling(nrow(z_sample)/batch_size))
        }
    }
    
    
    ### post-processing
    message("post-processing")
    # denormalization
    if (!is.null(minmax)) {
        x_gen <- gradDescent::minmaxDescaling(data.frame(x_gen), minmax$scalingParameter)
        x_train <- gradDescent::minmaxDescaling(data.frame(x_train), minmax$scalingParameter)
    }
    
    # original label
    if (!is.null(y_train)) {
        lenc <- x$preprocessing$lenc
        y_gen <- CatEncoders::inverse.transform(lenc, as.numeric(y_gen))
        y_train <- CatEncoders::inverse.transform(lenc, y_train)
    }
    
    result$x_gen <- x_gen
    result$y_gen <- y_gen
    result$x_train <- x_train
    result$y_train <- y_train
    result$latent_vector <- x_train_encoded
    result
}





plot_aug <- function(x, plot_fun, ...) {
    x_gen <- x$x_gen
    x_train <- x$x_train
    if (!is.null(x$y_gen)) {
        y_gen <- x$y_gen
        y_train <- x$y_train
    }
    
    # single cell object
    if (!is.null(x$y_gen)) {
        group <- data.frame(c(as.vector(y_train), as.vector(y_gen)),
                            c(rep("real", nrow(x_train)), rep("gen", nrow(x_gen))))
        rownames(group) <- c(rownames(x_train), paste("gen", seq_len(nrow(x_gen)), sep = ""))
        colnames(group) <- c("label", "group")
    } else {
        group <- data.frame(c(rep("real", nrow(x_train)), rep("gen", nrow(x_gen))))
        rownames(group) <- c(rownames(x_train), paste("gen", seq_len(nrow(x_gen)), sep = ""))
        colnames(group) <- "group"
    }
    
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(
            counts = t(as.matrix(rbind(as.matrix(x_train), as.matrix(x_gen))))
        ), 
        colData = group
    )
    sce <- scater::logNormCounts(sce)
    
    # plot
    if (!is.null(x$y_gen)) {
        do.call(eval(parse(text = paste("scater::plot", plot_fun, sep = ""))),
                list(do.call(eval(parse(text = paste("scater::run", plot_fun, sep = ""))),
                            list(sce)), colour_by = "group", shape_by = "label", ...))
    } else {
        do.call(eval(parse(text = paste("scater::plot", plot_fun, sep = ""))),
                list(do.call(eval(parse(text = paste("scater::run", plot_fun, sep = ""))),
                            list(sce)), colour_by = "group", ...))
    }
}