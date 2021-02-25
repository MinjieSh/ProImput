#' Line plot, compare all of the methods.
#' @param nrmse nrmse matrix for different methods and missing rates.
#' row: #repetition
#' column: #methods * #missing_rate
#' All methods as a group, one group per missing rate
#' @return line plot comparing imputation performance
#' of different methods at different missing rates.
plot_nrmse <-
  function(nrmse,
           method_name,
           missing_rate,
           ytitle,
           subtitle,
           panel,
           label = TRUE) {
    mean <- rowMeans(nrmse, na.rm = TRUE)
    sd <- apply(nrmse, 1, function(x)
      sd(x, na.rm = TRUE))
    
    number_of_missing_rate <- length(missing_rate)
    number_of_method <- length(method_name)
    
    methods = factor(method_name)
    methods = rep(methods, number_of_missing_rate)
    
    missingRate = missing_rate
    missingRate = rep(missingRate, each = number_of_method)
    
    df <- data.frame(missingRate, methods, mean, sd)
    
    res <- ggplot(df,
                  aes(
                    x = missingRate,
                    y = mean,
                    group = methods,
                    color = methods
                  )) +
      geom_line(size = 1, alpha = 0.8,
                aes(linetype = methods)) +
      geom_point(aes(shape = methods), size = 3) +
      scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 15, 16, 17, 18)) +
      xlab("Total Missing Rate") +
      ylab(ytitle) +
      ggtitle(subtitle) +
      coord_cartesian(xlim = c(
        missing_rate[1] - 1.5 * (missing_rate[2] - missing_rate[1]),
        missing_rate[number_of_missing_rate] + 1.5 * (missing_rate[2] -
                                                        missing_rate[1])
      )) +
      facet_grid(factor(
        mean < panel,
        level = c(FALSE, TRUE),
        labels = c("A", "B")
      ) ~ .,
      scales = "free_y")
    
    if (label) {
      return(
        res + geom_dl(aes(label = methods),
                      method = list(
                        dl.combine("first.points", "last.points"),
                        cex = 1
                      )) +
          theme(
            title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank()
          )
      )
    } else {
      return(
        res + theme(
          title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.text.x = element_text(color = "white"),
          axis.text.y = element_text(color = "white"),
          legend.text = element_text(size = 12, colour = "white"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_blank()
        )
      )
    }
  }

plot_missing_mechanism <- function(data, subtitle, label = TRUE) {
  missing_rate <- apply(data, 2, get_misg_rate)
  mean <- colMeans(data, na.rm = TRUE)
  
  mechanism <- ggplot(as.data.frame(missing_rate)) +
    geom_point(aes(x = mean, y = missing_rate), colour = "dodgerblue") +
    coord_cartesian(xlim = c(0, 25), ylim = c(0, 100)) +
    theme(legend.position = "none") +
    xlab("Mean log2 protein abundance") +
    ylab("Missing rate by protein (%)") +
    ggtitle(subtitle)
  
  if (label) {
    return(mechanism + theme(
      title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face =
                                  "bold")
    ))
  } else {
    return(
      mechanism + theme(
        title = element_text(
          size = 12,
          face = "bold",
          colour = "white"
        ),
        axis.title = element_text(
          size = 12,
          face = "bold",
          colour = "white"
        ),
        axis.text.x = element_text(color = "white"),
        axis.text.y = element_text(color = "white"),
        legend.text = element_text(size = 12, colour = "white")
      )
    )
  }
}

plot_pattern_change <-
  function(before, after, subtitle, label = TRUE) {
    before_missing_rate <- apply(before, 2, get_misg_rate)
    before_mean <- colMeans(before, na.rm = TRUE)
    
    after_missing_rate <- apply(after, 2, get_misg_rate)
    after_mean <- colMeans(after, na.rm = TRUE)
    
    group <- c(rep(1, ncol(before)), rep(2, ncol(after)))
    mean <- c(before_mean, after_mean)
    missing_rate <- c(before_missing_rate, after_missing_rate)
    df <- as.data.frame(cbind(group, mean, missing_rate))
    
    mechanism <- ggplot(df) +
      geom_point(aes(x = mean, y = missing_rate, color = group)) +
      coord_cartesian(xlim = c(0, 25), ylim = c(0, 100)) +
      theme(legend.position = "none") +
      xlab("Mean log2 protein abundance") +
      ylab("Missing rate by protein (%)") +
      ggtitle(subtitle)
    
    if (label) {
      return(mechanism + theme(
        title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face =
                                    "bold")
      ))
    } else{
      return(
        mechanism + theme(
          title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.text.x = element_text(color = "white"),
          axis.text.y = element_text(color = "white"),
          legend.text = element_text(size = 12, colour = "white")
        )
      )
    }
    
  }

plot_SOR <-
  function(sor_matrix,
           method_name,
           missing_rate,
           subtitle,
           label = TRUE) {
    sor <- matrix(log10(t(sor_matrix)),
                  nrow = ncol(sor_matrix) * nrow(sor_matrix),
                  ncol = 1)
    
    number_of_missing_rate <- length(missing_rate)
    number_of_method <- length(method_name)
    
    methods = factor(method_name)
    methods = rep(methods, number_of_missing_rate)
    
    missingRate = missing_rate
    missingRate = rep(missingRate, each = number_of_method)
    
    df <- data.frame(missingRate, methods, sor)
    
    res <- ggplot(df,
                  aes(
                    x = missingRate,
                    y = sor,
                    group = methods,
                    color = methods
                  )) +
      geom_line(size = 1, alpha = 0.8,
                aes(linetype = methods)) +
      geom_point(aes(shape = methods), size = 3) +
      scale_shape_manual(values = c(0, 1, 2, 15, 17, 5, 6, 16)) +
      xlab("Total Missing Rate") +
      ylab("log10 NRMSE-based Sum of Ranks (SOR)") +
      ggtitle(subtitle) +
      coord_cartesian(xlim = c(
        missing_rate[1] - 1.5 * (missing_rate[2] - missing_rate[1]),
        missing_rate[number_of_missing_rate] + 1.5 * (missing_rate[2] -
                                                        missing_rate[1])
      ))
    
    
    if (label) {
      return(
        res + geom_dl(aes(label = methods),
                      method = list(
                        dl.combine("first.points", "last.points"),
                        cex = 1
                      )) +
          theme(
            title = element_text(size = 12, face = "bold"),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12)
          )
      )
    } else{
      return(
        res + theme(
          title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.title = element_text(
            size = 12,
            face = "bold",
            colour = "white"
          ),
          axis.text.x = element_text(color = "white"),
          axis.text.y = element_text(color = "white"),
          legend.text = element_text(size = 12, colour = "white")
        )
      )
    }
  }
