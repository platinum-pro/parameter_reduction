library(ggplot2)
library(dplyr)
library(gridExtra)
library(purrr)
library(tidyr)

# Function to calculate RMSE
calculate_rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

# Create the data frame - only OpenH
#data_1e <- data.frame(
 # x = c(0, 0.05, 0.13, 0.25, 0.5, 1, 3, 6, 11, 35, 70, 140, 280, 560, 1120),
  #OpenH = c(95.37, 95.32, 95.32, 95.32, 95.32, 94.735, 94.725, 94.725, 91.215, 73.065, 44.385, 22.8, 5.845, 2.93, 2.92)
#)

#conditions <- c('OpenH')
#condition_names <- c('Aggregate Demand for Health Insurance')

# Create the data frame - only OpenH
data_1e <- data.frame(
  x = c(0, 10, 25, 50, 75, 100, 150, 250, 400, 550, 700, 850, 1000, 1500, 2000),
  OpenH = c(82.774365, 77.283636, 71.772125, 62.973932, 55.091854, 48.742754, 42.845601, 36.430512, 27.437269, 22.465737, 18.759734, 18.426704, 15.445512, 13.943189, 13.671014)
)

conditions <- c('OpenH')
condition_names <- c('Aggregate Demand for HIV Vaccines')


# Model functions
KNK <- function(x, alpha, Qo) {
  Qo * 10^(exp(-alpha * Qo * x) - 1)
}

KNKEQ <- function(x, alpha, condition) {
  Qo <- as.numeric(data_1e[data_1e$x == 0, condition])
  if(length(Qo) == 0 || !is.numeric(Qo)) {
    stop("Invalid Qo value for condition: ", condition)
  }
  Qo * 10^(exp(-alpha * Qo * x) - 1)
}

# Modified Koff function with fixed k
Koff <- function(x, alpha, Qo) {
  k <- log10(max(data_1e$OpenH)) - log10(min(data_1e$OpenH))
  Qo * 10^(k * exp(-alpha * Qo * x) - 1)
}

# Updated fit_model function
fit_model <- function(data, condition) {
  Qo_initial <- as.numeric(data[data$x == 0, condition])
  if(length(Qo_initial) == 0 || !is.numeric(Qo_initial)) {
    stop("Invalid initial Qo value for condition: ", condition)
  }
  
  fits <- list()
  
  # KNK
  formula1 <- as.formula(paste(condition, "~ KNK(x, alpha, Qo)"))
  fits[["KNK"]] <- tryCatch({
    nls(formula = formula1, data = data, start = list(alpha = 0.000001, Qo = Qo_initial),
        algorithm = "port", lower = c(alpha = 0, Qo = 0), upper = c(alpha = Inf, Qo = Inf))
  }, error = function(e) {
    message("Error fitting KNK model: ", e$message)
    NULL
  })
  
  # KNKEQ
  formula2 <- as.formula(paste(condition, "~ KNKEQ(x, alpha, '", condition, "')", sep = ""))
  fits[["KNKEQ"]] <- tryCatch({
    nls(formula = formula2, data = data, start = list(alpha = 0.000001),
        algorithm = "port", lower = c(alpha = 0), upper = c(alpha = Inf))
  }, error = function(e) {
    message("Error fitting KNKEQ model: ", e$message)
    NULL
  })
  
  # Modified Koff
  formula3 <- as.formula(paste(condition, "~ Koff(x, alpha, Qo)"))
  fits[["Koff"]] <- tryCatch({
    nls(formula = formula3, 
        data = data, 
        start = list(alpha = 0.000001, Qo = 100),
        algorithm = "port",
        lower = c(alpha = 0, Qo = 0), 
        upper = c(alpha = Inf, Qo = Inf),
        control = nls.control(maxiter = 50000))
  }, error = function(e) {
    message("Error fitting Koff model: ", e$message)
    NULL
  })
  
  fits
}

# Bootstrapping function
bootstrap_analysis <- function(data, condition, fit, model_name, n_iterations = 1000) {
  if (is.null(fit)) return(list(results = NULL, successful_iterations = 0))
  
  fit_formula <- formula(fit)
  
  results <- matrix(NA, nrow = n_iterations, ncol = length(coef(fit)), 
                    dimnames = list(NULL, names(coef(fit))))
  
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    tryCatch({
      bootstrap_sample <- data[sample(nrow(data), replace = TRUE), ]
      bootstrap_fit <- nls(fit_formula, data = bootstrap_sample, 
                           start = coef(fit),
                           algorithm = "port",
                           lower = rep(0, length(coef(fit))),
                           upper = rep(Inf, length(coef(fit))))
      results[i, ] <- coef(bootstrap_fit)
      successful_iterations <- successful_iterations + 1
    }, error = function(e) {
      message(sprintf("Error in iteration %d of bootstrap analysis for condition '%s', model '%s': %s", 
                      i, condition, model_name, e$message))
    })
  }
  
  list(results = as.data.frame(na.omit(results)), successful_iterations = successful_iterations)
}

# Create plot function
create_plot <- function(condition, fits, data, condition_name) {
  if (all(sapply(fits, is.null))) return(NULL)
  
  x_range <- 10^seq(log10(0.01), log10(max(data$x)), length.out = 100)
  
  plot_data <- data.frame(x = x_range)
  model_names <- c("KNK", "KNKEQ", "Koff")
  colors <- c("black", "black", "black")
  line_types <- c("dotted", "dashed", "solid")
  
  for (model_name in model_names) {
    if (!is.null(fits[[model_name]])) {
      plot_data[[model_name]] <- predict(fits[[model_name]], newdata = list(x = x_range))
    }
  }
  
  plot_data_long <- tidyr::pivot_longer(plot_data, 
                                        cols = intersect(names(plot_data), model_names), 
                                        names_to = "Model", 
                                        values_to = "Prediction")
  
  p <- ggplot() +
    geom_point(data = data, aes_string(x = "x", y = condition)) +
    geom_line(data = plot_data_long, aes(x = x, y = Prediction, color = Model, linetype = Model)) +
    scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100, 1000), labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
    labs(title = condition_name, x = "Price ($)", y = "Average likelihood \n of accepting HIV vaccine (%)") +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.line = element_line(color = "black"),
          plot.title = element_text(hjust = 0.5)) +
    annotation_logticks(sides = "b") +
    scale_color_manual(values = colors[1:length(unique(plot_data_long$Model))],
                       labels = gsub("_", " ", unique(plot_data_long$Model))) +
    scale_linetype_manual(values = line_types[1:length(unique(plot_data_long$Model))],
                          labels = gsub("_", " ", unique(plot_data_long$Model)))
  
  p
}

# Function to display results
display_results <- function(condition, fits, bootstrap_results, condition_name) {
  cat("\nResults for", condition_name, ":\n")
  
  model_names <- c("KNK", "KNKEQ", "Koff")
  
  for (model_name in model_names) {
    cat("\n", model_name, ":\n")
    if (!is.null(fits[[model_name]])) {
      params <- coef(fits[[model_name]])
      cat("Parameters:\n")
      print(params)
      
      # Calculate R-squared
      residuals <- residuals(fits[[model_name]])
      tss <- sum((data_1e[[condition]] - mean(data_1e[[condition]], na.rm = TRUE))^2, na.rm = TRUE)
      rss <- sum(residuals^2)
      r_squared <- 1 - (rss / tss)
      cat("R-squared =", round(r_squared, 4), "\n")
      
      # Calculate RMSE
      predicted_values <- predict(fits[[model_name]])
      rmse <- sqrt(mean((data_1e[[condition]] - predicted_values)^2))
      cat("RMSE =", round(rmse, 4), "\n")
      
      # Calculate AIC and BIC
      aic <- AIC(fits[[model_name]])
      bic <- BIC(fits[[model_name]])
      cat("AIC =", round(aic, 4), "\n")
      cat("BIC =", round(bic, 4), "\n")
      
      # Display bootstrapping results
      if (!is.null(bootstrap_results[[model_name]]$results)) {
        cat("\nBootstrapping Analysis Summary:\n")
        print(summary(bootstrap_results[[model_name]]$results))
        
        # Calculate confidence intervals
        ci <- apply(bootstrap_results[[model_name]]$results, 2, quantile, probs = c(0.025, 0.975))
        cat("\n95% Confidence Intervals:\n")
        print(ci)
        
        cat("\nSuccessful bootstrap iterations:", bootstrap_results[[model_name]]$successful_iterations, "\n")
      }
    } else {
      cat("Model fitting failed for this condition.\n")
    }
  }
}

# Apply analysis to all conditions
results <- map(conditions, ~fit_model(data_1e, .x))

# Perform bootstrapping
bootstrap_results <- map2(conditions, results, function(condition, fits) {
  map2(fits, names(fits), ~bootstrap_analysis(data_1e, condition, .x, .y))
})

# Create plots
plots <- map2(conditions, results, ~create_plot(.x, .y, data_1e, condition_names[which(conditions == .x)]))

# Arrange plots in a grid
grid.arrange(grobs = plots[!sapply(plots, is.null)], ncol = 1)

# Display results for all conditions
walk2(conditions, results, ~display_results(.x, .y, 
                                            bootstrap_results[[which(conditions == .x)]], 
                                            condition_names[which(conditions == .x)]))

# Calculate and display total successful bootstrap iterations
total_successful_iterations <- sum(unlist(map(bootstrap_results, ~map_dbl(.x, ~.x$successful_iterations))))
cat("\nTotal successful bootstrap iterations across all conditions and models:", total_successful_iterations, "\n")