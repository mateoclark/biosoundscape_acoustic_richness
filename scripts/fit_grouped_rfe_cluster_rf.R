fit_grouped_rfe_cluster_rf <- function(cor_data, predictor_vars, metric,
                                       cluster_col = "CLUSTER_ID",
                                       recnum_col = "LogRecNum",
                                       top_n_per_group = 3,
                                       rfe_sizes = seq(10, 100, by = 10),
                                       seed = 42,
                                       num.trees = 500,
                                       min.node.size = 5,
                                       num.threads = 2,
                                       rfe_on = 1,
                                       importance = "impurity",
                                       verbose = TRUE) {
  
  set.seed(seed)
  
  # Step 1: Prepare data
  cor_data <- cor_data %>% rename(response = !!sym(metric))
  df_rf <- cor_data %>%
    dplyr::select(response, all_of(predictor_vars), all_of(cluster_col), all_of(recnum_col)) %>%
    filter(!is.na(response))
  df_rf[is.na(df_rf)] <- 0
  df_rf <- df_rf %>%
    mutate(cluster_id = .data[[cluster_col]],
           LogRecNum = .data[[recnum_col]]) %>%
    select(-all_of(cluster_col))
  
  group_labels <- sub("^[^_]+_", "", predictor_vars)
  names(group_labels) <- predictor_vars

  # Step 2: Full model for variable importance
  full_model <- ranger(
    response ~ .,
    data = df_rf %>% select(-cluster_id,-LogRecNum),
    importance = importance,
    mtry = floor(sqrt(ncol(df_rf) - 2)),
    num.trees = num.trees,
    min.node.size = min.node.size,
    num.threads = num.threads,
    seed = seed
  )
  
  imp <- importance(full_model)
  imp_df <- data.frame(
    Variable = names(imp),
    Importance = as.numeric(imp),
    Group = group_labels[names(imp)]
  )
  
  # Step 3: Group-wise top-N filtering
  top_vars <- imp_df %>%
    group_by(Group) %>%
    slice_max(order_by = Importance, n = top_n_per_group, with_ties = FALSE) %>%
    pull(Variable)
  
  filtered_data <- df_rf %>% dplyr::select(all_of(top_vars), response, cluster_id)
  
  # Step 4: Drop highly correlated variables
  cor_mat <- cor(filtered_data %>% dplyr::select(-response, -cluster_id), use = "pairwise.complete.obs")
  high_corr <- findCorrelation(cor_mat, cutoff = 0.9, names = TRUE)
  final_vars <- setdiff(colnames(filtered_data)[!colnames(filtered_data) %in% c("response", "cluster_id")], high_corr)
  
  # Step 5: Prepare final data and define custom CV folds
  df_final <- filtered_data %>% dplyr::select(all_of(final_vars), response, cluster_id)
  
  # Convert CLUSTER_ID to character
  df_final$cluster_id <- as.character(df_final$cluster_id)
  df_final$cluster_id[df_final$cluster_id == "-1"] <- "Unclustered"
  
  # Set up cluster folds
  group_folds <- groupKFold(df_final$cluster_id, k = 10)
  train_index <- group_folds
  
  # Drop cluster column before modeling
  df_final_model <- df_final %>% select(-cluster_id)
  
  # Step 6: RFE (if enabled)
  if (rfe_on == 1) {
    ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 3, allowParallel = TRUE)
    ctrl$functions$fit <- function(x, y, first, last, ...) {
      train(x, y, method = "ranger", importance = importance, metric = "Rsquared", num.threads = 1, ...)
    }
    
    rfe_result <- rfe(
      x = df_final_model %>% select(-response),
      y = df_final_model$response,
      sizes = rfe_sizes[rfe_sizes <= length(final_vars)],
      rfeControl = ctrl
    )
    
    selected_vars <- rfe_result$optVariables
  } else {
    selected_vars <- final_vars
    rfe_result <- NULL
  }
  
  # Step 7: Prepare final modeling data
  final_data <- df_final_model %>% dplyr::select(all_of(selected_vars), response)
  if (metric == "AIrichness") final_data$LogRecNum <- df_rf$LogRecNum
  
  # Step 8: Tune mtry
  tune_ctrl <- trainControl(
    method = "cv",
    number = 10,
    index = train_index,
    allowParallel = TRUE
  )
  
  tuned_model <- train(
    response ~ .,
    data = final_data,
    method = "ranger",
    trControl = tune_ctrl,
    tuneGrid = expand.grid(
      mtry = seq(1, length(selected_vars), by = 1),
      splitrule = "variance",
      min.node.size = min.node.size
    ),
    importance = importance,
    metric = "Rsquared",
    num.threads = 1
  )
  
  best_mtry <- tuned_model$bestTune$mtry
  fold_results <- tuned_model$resample
  
  # Step 9: Final model
  final_model <- ranger(
    response ~ .,
    data = final_data,
    mtry = best_mtry,
    splitrule = "variance",
    min.node.size = min.node.size,
    num.trees = num.trees,
    importance = importance,
    num.threads = num.threads,
    seed = seed + 3
  )
  
  if (verbose) {
    cat("Selected", length(selected_vars), "variables\n")
    cat("Final model R\u00b2:", round(final_model$r.squared, 4), "\n")
  }
  
  list(
    full_model    = full_model,
    final_model   = final_model,
    selected_vars = selected_vars,
    rfe_result    = rfe_result,
    final_data    = final_data,
    fold_results  = fold_results,
    metric        = metric
  )
}
