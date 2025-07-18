#' Selecting genotypes in multiple environments
#' 
#' @param data         A data.frame containing "Taxa", "Env", and "Obs".
#' @param K            An additive genetic relationship matrix.
#' @param n            Size of the recommendation candidate set.
#' @param method       "maximum", "minimum" or "target_value"
#' @param target_value Target value used when method is "target_value".
#' @export
#' @import sommer
#' @import MASS
#' @import dplyr
#' @import ggplot2
#' @examples
#'

#' del <- sample(1359,500)
#' data <- Trait
#' data$Obs[del] <- NA
#' result <- select_genotypes(data, K, n=c(20,30,40), method = "maximum")



select_genotypes <- function(data, K, n, method = c("maximum", "minimum", "target_value"), target_value = NULL) {

  k <- nrow(K)
  data <- data %>%
    arrange(Env, Taxa)
  
  
  # Fit the MGE model
  MGE <- tryCatch({
    mmes(
      Obs ~ Env,
      random = ~ vsm(ism(Taxa), Gu = K) + vsm(dsm(Env), ism(Taxa), Gu = K),
      rcov = ~ vsm(dsm(Env), ism(units)),
      data = data,
      verbose = FALSE,
      tolParInv = 1e-06
    )
  }, error = function(e) {
    for (rep in 1:50) {
      tolParInv <- 1e-06 + 0.005 * rep
      MGE <- tryCatch({
        mmes(
          Obs ~ Env,
          random = ~ vsm(ism(Taxa), Gu = K) + vsm(dsm(Env), ism(Taxa), Gu = K),
          rcov = ~ vsm(dsm(Env), ism(units)),
          data = data,
          verbose = FALSE,
          tolParInv = tolParInv
        )
      }, error = function(e) NULL)
      if (!is.null(MGE)) break
    }
    MGE
  })

  if (is.null(MGE)) {
    message("MGE model is NULL. Computation stopped.")
    return()
  }
  
  # Predict GEBVs
  Z0 <- rep(c(MGE[["b"]][1], MGE[["b"]][1] + MGE[["b"]][-1]), each = nrow(K))
  Z1 <- unlist(MGE[["u"]][1:k]) + unlist(MGE[["u"]][-(1:k)])
  result <- data.frame(GEBV = (Z0 + Z1), Taxa = rep(colnames(K), 3), Env = rep(1:3, each = nrow(K)))
  
  # Select genotypes
  output <- data.frame()
  for (i in seq_along(n)) {
    sub_result <- result[data$Env == i & is.na(data$Obs), ]
    
    if (method == "maximum") {
      ordered <- sub_result[order(-sub_result$GEBV), ]
    } else if (method == "minimum") {
      ordered <- sub_result[order(sub_result$GEBV), ]
    } else if (method == "target_value") {
      if (is.null(target_value)) stop("Please specify a target_value when using method = 'target_value'")
      ordered <- sub_result[order(abs(sub_result$GEBV - target_value)), ]
    }
    
    selected <- ordered[1:n[i], ]
    output <- rbind(output, selected)
  }
  final <- list()
  final[['recommend']] <- output
  
  # Plot  
  plot_tbv_lines <- function(data, n_group = 8) {
    
    df <- matrix(result$GEBV,ncol=max(result$Env)) 
    colnames(df) <- paste0("ENV", seq_len(max(result$Env)))
    df <- data.frame(df)
    
    df$line_id <- factor(1:k)
    
    df$group <- cut(df$ENV1, breaks = quantile(df$ENV1, probs = seq(0, 1, length.out = n_group + 1)), 
                    include.lowest = TRUE, labels = paste0("G", 1:n_group))
    
    df_long <- tidyr::pivot_longer(df, cols = starts_with("ENV"),
                                   names_to = "Environment", values_to = "Value")
    df_long$Environment <- factor(df_long$Environment,
                                  levels = paste0("ENV", seq_len(max(result$Env))))
    
    p <- ggplot(df_long, aes(x = Environment, y = Value, group = line_id, color = group)) +
      geom_line(size = 0.8, alpha = 0.7) +
      geom_point(size = 1.2, alpha = 0.5) +
      labs(
        title = 'G×E plot',
        y = "GEBV",
        x = "Environment"
      ) +
      theme_minimal(base_family = "serif") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none"
      )
    
    return(p)
  }
  final[['plot']] <- plot_tbv_lines(result)
  return(final)
}
