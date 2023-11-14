#' @export
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices devAskNewPage
#' @importFrom grDevices palette
#' @importFrom stats aggregate
plot.rarecol <- function(x, xlim1, xlim2, ylim, ylim1, ylim2, palette1,
                         palette2, main1, main2, xlab1, xlab2, ylab1, ylab2,
                         las1 = 1, las2 = 1, cextText = 0.75,
                         legendbar = TRUE, ...) {

    # Prepare genetic curves data
    data_gen <- x[[1]]
    n_pop <- sort(unique(data_gen$Populations))
    n_pos <- sort(unique(data_gen$VariablePositions))
    data_sum_gen <- aggregate(x   = data_gen,
                              by  = list(data_gen$GeneticReplicate,
                                         data_gen$VariablePositions,
                                         data_gen$Populations),
                              FUN = mean)
    data_sum_gen <- data_sum_gen[, -c(1:3)]
    data_sum_gen <- aggregate(x   = data_sum_gen,
                              by  = list(data_sum_gen$VariablePositions,
                                         data_sum_gen$Populations),
                              FUN = mean)
    data_sum_gen <- data_sum_gen[, -c(1:2)]

    # Prepare field curves data
    data_field <- x[[2]]
    data_sum_field <- aggregate(x   = data_field,
                                by  = list(data_field$FieldReplicate,
                                           data_field$VariablePositions,
                                           data_field$Populations),
                                FUN = mean)
    data_sum_field <- data_sum_field[, -c(1:3)]
    data_sum_field <- aggregate(x   = data_sum_field,
                                by  = list(data_sum_field$VariablePositions,
                                           data_sum_field$Populations),
                                FUN = mean)
    data_sum_field <- data_sum_field[, -c(1:2)]

    # Predetermined arguments
    if (missing(palette1)) {
        palette <- colorRampPalette(c("blue", "green", "yellow", "red"))
        palette1 <- palette(length(n_pos))
    }
    if (missing(palette2)) {
        palette <- colorRampPalette(c("blue", "green", "yellow", "red"))
        palette2 <- palette(length(n_pop))
    }
    if (missing(xlim1)) {
        if(legendbar == TRUE) {
            max <- max(data_sum_gen$Populations) * 1.3
        } else {
            max <- max(data_sum_gen$Populations)
        }
        min <- 0
        xlim1 <- c(min, max)
    }
    if (missing(xlim2)) {
        if (legendbar == TRUE) {
            max <- max(data_sum_field$VariablePositions) * 1.3
        } else {
            max <- max(data_sum_field$VariablePositions)
        }
        min <- 0
        xlim2 <- c(min, max)
    }
    if (missing(ylim1)) {
        if (missing(ylim)) {
            min <- 0
            max <- max(data_sum_gen$TotalColonizationEvents,
                       data_sum_field$TotalColonizationEvents)
            ylim1 <- c(min, max)
        } else {
            ylim1 <- ylim
        }
    }
    if (missing(ylim2)) {
        if (missing(ylim)) {
            min <- 0
            max <- max(data_sum_gen$TotalColonizationEvents,
                       data_sum_field$TotalColonizationEvents)
            ylim2 <- c(min, max)
        } else {
            ylim2 <- ylim
        }
    }
    if (missing(main1)) {
        main1 <- "Genetic estimator"
    }
    if (missing(xlab1)) {
        if (length(unique(data_field[[1]])) ==
            length(unique(data_field[[2]]))) {
            
            xlab1 <- "Individuals"
        } else {
            xlab1 <- "Populations"
        }
    }
    if (missing(ylab1)) {
        ylab1 <- "Colonization events"
    }
    if (missing(main2)) {
        main2 <- "Field estimator"
    }
    if (missing(xlab2)) {
        xlab2 <- "Variable positions"
    }
    if (missing(ylab2)) {
        ylab2 <- "Colonization events"
    }

    # Plot 1
    px <- data_sum_gen$Populations[data_sum_gen$VariablePositions == n_pos[1]]
    py <- data_sum_gen$TotalColonizationEvents[data_sum_gen$VariablePositions ==
                                               n_pos[1]]
    plot(x    = px,
         y    = py,
         xlim = xlim1,
         ylim = ylim1,
         col  = palette1[1],
         main = main1,
         xlab = xlab1,
         ylab = ylab1,
         type = "l",
         las  = las1, ...)
    if (legendbar == TRUE) {
        box_min <- (ylim1[2] - ylim1[1]) * 0.2 + ylim1[1]
        box_max <- (ylim1[2] - ylim1[1]) * 0.8 + ylim1[1]
        text(x      = xlim1[2] * 0.85,
             y      = (ylim1[2] - ylim1[1]) / 2 + ylim1[1],
             labels = "Variable positions",
             cex    = cextText,
             srt    = 90)
        text(x      = xlim1[2],
             y      = (ylim1[2] - ylim1[1]) * 0.1 + ylim1[1],
             labels = min(n_pos),
             adj    = 1,
             cex    = cextText)
        text(x      = xlim1[2],
             y      = (ylim1[2] - ylim1[1]) * 0.9 + ylim1[1],
             labels = max(n_pos),
             adj    = 1,
             cex    = cextText)
    }
    for (i in 1:length(n_pos)) {
        if (i > min(n_pos)) {
            lx <- data_sum_gen$Populations
            ly <- data_sum_gen$TotalColonizationEvents
            lines(x   = lx[data_sum_gen$VariablePositions == n_pos[i]],
                  y   = ly[data_sum_gen$VariablePositions == n_pos[i]],
                  col = palette1[i])
        }
        if (legendbar == TRUE) {
            rect(xleft   = xlim1[2] * 0.92,
                 ybottom = box_min + (box_max - box_min)/length(n_pos) *
                           (length(n_pos[1:i]) - 1),
                 xright  = xlim1[2],
                 ytop    = box_min + (box_max - box_min)/length(n_pos) *
                           length(n_pos[1:i]),
                 col     = palette1[i],
                 border  = NA)
        }
    }

    if (prod(par("mfrow")) <= 2) devAskNewPage(ask = TRUE)

    # Plot 2
    plot2x <- data_sum_field$VariablePositions
    plot2y <- data_sum_field$TotalColonizationEvents
    plot(x    = plot2x[data_sum_field$Populations == n_pop[1]],
         y    = plot2y[data_sum_field$Populations == n_pop[1]],
         main = main2,
         xlab = xlab2,
         ylab = ylab2,
         type = "l",
         col  = palette2[1],
         ylim = ylim2,
         xlim = xlim2,
         las  = las2, ...)
    if (legendbar == TRUE) {
        box_min <- (ylim2[2] - ylim2[1]) * 0.2 + ylim2[1]
        box_max <- (ylim2[2] - ylim2[1]) * 0.8 + ylim2[1]
        text(x      = xlim2[2] * 0.85,
             y      = (ylim2[2] - ylim2[1]) / 2 + ylim2[1],
             labels = xlab1,
             cex    = cextText,
             srt    = 90)
        text(x      = xlim2[2],
             y      = (ylim2[2] - ylim2[1]) * 0.1 + ylim2[1],
             labels = min(n_pop),
             adj    = 1,
             cex    = cextText)
        text(x      = xlim2[2],
             y      = (ylim2[2] - ylim2[1]) * 0.9 + ylim2[1],
             labels = max(n_pop),
             adj    = 1,
             cex    = cextText)
    }
    for (i in n_pop) {
        if (i > min(n_pop)) {
            lx <- data_sum_field$VariablePositions
            ly <- data_sum_field$TotalColonizationEvents
            lines(x   = lx[data_sum_field$Populations == i],
                  y   = ly[data_sum_field$Populations == i],
                  col = palette2[i])
        }
        if (legendbar == TRUE) {
            rect(xleft   = xlim2[2] * 0.92,
                 ybottom = box_min + (box_max - box_min)/length(n_pop) *
                           (length(n_pop[1:i]) - 1),
                 xright  = xlim2[2],
                 ytop    = box_min + (box_max - box_min)/length(n_pop) *
                           length(n_pop[1:i]),
                 col     = palette2[i],
                 border  = NA)
        }
    }
}
