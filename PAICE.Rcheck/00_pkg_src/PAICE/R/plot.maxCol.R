#' @export
#' @importFrom graphics arrows
#' @importFrom graphics axis
#' @importFrom graphics box
#' @importFrom graphics curve
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics plot.default
#' @importFrom graphics points
#' @importFrom graphics rect
#' @importFrom graphics title
plot.maxCol <- function(x, xlim, ylim, col, xlabbotton, xlabtop, ylab,
                        main, pch = 16, lty = 1, lwd = 2, cex = 1,
                        estimation = TRUE, legend = TRUE, ...) {
    if (is.numeric(x$ParametersGen)) {
        f_gen <- function(x) x$ParametersGen[1, 1] * x /
                             (x$ParametersGen[2, 1] + x) +
                             x$ParametersGen[3, 1]
    } else {
        f_gen <- NULL
    }
    if (is.numeric(x$ParametersField)) {
        f_field <- function(x) x$ParametersField[1, 1] * (x - 1) /
                               (x$ParametersField[2, 1] + x - 1)
    } else {
        f_field <- NULL
    }
    if (missing(main)) {
        main = "Estimated number of colonization events"
    }
    if (missing(xlabbotton)) {
        xlabbotton <- "Variable positions"
    }
    if (missing(xlabtop)) {
        if (names(x$DataField)[2] == "populations") {
            xlabtop <- "Populations"
        } else {
            xlabtop <- "Individuals"
        }
    }
    if (missing(ylab)) {
        ylab <- "Estimated colonization events"
    }
    if (missing(xlim)) {
        if (estimation == TRUE) {
            xlim <- c(0, max(x$DataGen$positions,
                             x$DataField[[2]] * 1.2))
        } else {
            xlim <- c(0, max(x$DataGen$positions, x$DataField[[2]]))
        }
    }
    if (missing(ylim)) {
        if (estimation == TRUE & is.numeric(x$Summary)) {
            ylim <- c(0, max(x$Summary,
                             na.rm = TRUE) * 1.1)
        } else {
            ylim <- c(0, max(x$DataGen$colonization,
                             x$DataField$colonization,
                             na.rm = TRUE))
        }
    }

    if (missing(col)) {
        col <- c("#f7746a", "#36acae")
    }

    # Set margin settings
    original_mar <- par()
    new_mar <- original_mar$mar + c(0, 0, 1, 0)
    par(mar = new_mar)
    
    plot.default(x    = 1,
                 type = "n",
                 xlim = xlim,
                 ylim = ylim,
                 xlab = "",
                 ylab = ylab,
                 axes = FALSE, ...)
    axis(side      = 1,
         pos       = ylim[1] - (ylim[2] - ylim[1]) * 0.04,
         labels    = TRUE,
         las       = 1,
         at        = round(seq(from       = min(x$DataGen[[2]]),
                               to         = max(x$DataGen[[2]]),
                               length.out = 4)),
         col.axis  = col[1],
         col.ticks = col[1])
    axis(side   = 2,
         pos    = xlim[1] - (xlim[2] - xlim[1]) * 0.04,
         labels = TRUE,
         las    = 1)
    axis(side      = 3,
         pos       = ylim[2] + (ylim[2] - ylim[1]) * 0.04,
         labels    = TRUE,
         las       = 1,
         at        = round(seq(from       = min(x$DataField[[2]]),
                               to         = max(x$DataField[[2]]),
                               length.out = 4)),
         col.axis  = col[2],
         col.ticks = col[2])
    rect(xleft   = xlim[1] - (xlim[2] - xlim[1]) * 0.04,
         ybottom = ylim[1] - (ylim[2] - ylim[1]) * 0.04,
         xright  = xlim[2] + (xlim[2] - xlim[1]) * 0.04,
         ytop    = ylim[2] + (ylim[2] - ylim[1]) * 0.04,
         col     = NULL,
         border  = "black")
    title(main = main,
          line = 4)
    box()
    if (estimation == TRUE) {
        par(xpd = TRUE)
        rect(xleft   = xlim[2] * 0.8958,
             ybottom = ylim[1] - ylim[2] * 0.1,
             xright  = xlim[2] * 0.9167,
             ytop    = ylim[2] * 1.1,
             col     = "white",
             border  = "white")
        lines(x = c(xlim[2] * 0.8958 - (xlim[2] - xlim[1]) * 0.01,
                    xlim[2] * 0.8958 + (xlim[2] - xlim[1]) * 0.01),
              y = c(ylim[2] + (ylim[2] - ylim[1]) * 0.04 - (ylim[2] - ylim[1]) *
                    0.04,
                    ylim[2] + (ylim[2] - ylim[1]) * 0.04 + (ylim[2] - ylim[1]) *
                    0.04))
        lines(x = c(xlim[2] * 0.9167 - (xlim[2] - xlim[1]) * 0.01,
                    xlim[2] * 0.9167 + (xlim[2] - xlim[1]) * 0.01),
              y = c(ylim[2] + (ylim[2] - ylim[1]) * 0.04 - (ylim[2] - ylim[1]) *
                    0.04,
                    ylim[2] + (ylim[2] - ylim[1]) * 0.04 + (ylim[2] - ylim[1]) *
                    0.04))
        lines(x = c(xlim[2] * 0.8958 - (xlim[2] - xlim[1]) * 0.01,
                  xlim[2] * 0.8958 + (xlim[2] - xlim[1]) * 0.01),
             y = c(ylim[1] - (ylim[2] - ylim[1]) * 0.04 - (ylim[2] - ylim[1]) * 0.04,
                   ylim[1] - (ylim[2] - ylim[1]) * 0.04 + (ylim[2] - ylim[1]) * 0.04))
        lines(x = c(xlim[2] * 0.9167 - (xlim[2] - xlim[1]) * 0.01,
                  xlim[2] * 0.9167 + (xlim[2] - xlim[1]) * 0.01),
             y = c(ylim[1] - (ylim[2] - ylim[1]) * 0.04 - (ylim[2] - ylim[1]) * 0.04,
                   ylim[1] - (ylim[2] - ylim[1]) * 0.04 + (ylim[2] - ylim[1]) * 0.04))
        par(xpd = FALSE)
    }
    mtext(text = xlabbotton,
         side = 1,
         line = 3,
         col  = col[1])
    mtext(text = xlabtop,
         side = 3,
         line = 2.5,
         col  = col[2])
    points(x   = x$DataGen[[2]],
          y   = x$DataGen[[1]],
          pch = pch,
          col = col[1],
          cex = cex)
    if(is.null(f_gen) == FALSE) {
            curve(f_gen,
             from = 0,
             to   = max(x$DataGen[[2]]),
             col  = col[1],
             lty  = lty,
             lwd  = lwd,
             add  = TRUE)
    }
    points(x   = x$DataField[[2]],
           y   = x$DataField[[1]],
           pch = pch,
           col = col[2],
           cex = cex)
    if(is.null(f_field) == FALSE) {
        curve(f_field,
             from = 1,
             to   = max(x$DataField[[2]]),
             col  = col[2],
             lty  = lty,
             lwd  = lwd,
             add  = TRUE)
    }


    if(estimation == TRUE & is.numeric(x$Summary)) {
        points(x   = xlim[2] * 0.9583,
              y   = x$Summary[1, 1],
              col = col[1],
              pch = 16,
              cex = cex)
        if(is.na(x$Summary[1, 2])) {
            arrows(x0     = xlim[2] * 0.9583,
                  y0     = ylim[1],
                  x1     = xlim[2] * 0.9583,
                  y1     = x$Summary[1, 1],
                  length = 0.1,
                  angle  = 30,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[1])
        } else {
            arrows(x0     = xlim[2] * 0.9583,
                  y0     = x$Summary[1, 2],
                  x1     = xlim[2] * 0.9583,
                  y1     = x$Summary[1, 1],
                  length = 0.1,
                  angle  = 90,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[1])
        }
        if(is.na(x$Summary[1, 3])) {
            arrows(x0     = xlim[2] * 0.9583,
              y0     = ylim[2],
              x1     = xlim[2] * 0.9583,
              y1     = x$Summary[1, 1],
              length = 0.1,
              angle  = 30,
              lwd    = lwd,
              lty    = 1,
              code   = 1,
              col    = col[1])
        } else {
            arrows(x0     = xlim[2] * 0.9583,
                  y0     = x$Summary[1, 3],
                  x1     = xlim[2] * 0.9583,
                  y1     = x$Summary[1, 1],
                  length = 0.1,
                  angle  = 90,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[1])
        }
        points(x   = xlim[2],
              y   = x$Summary[2, 1],
              col = col[2],
              pch = 16,
              cex = cex)
        if(is.na(x$Summary[2, 2])) {
            arrows(x0     = xlim[2],
                  y0     = ylim[1],
                  x1     = xlim[2],
                  y1     = x$Summary[2, 1],
                  length = 0.1,
                  angle  = 30,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[2])
        } else {
            arrows(x0     = xlim[2],
                  y0     = x$Summary[2, 2],
                  x1     = xlim[2],
                  y1     = x$Summary[2, 1],
                  length = 0.1,
                  angle  = 90,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[2])
        }
        if(is.na(x$Summary[2, 3])) {
            arrows(x0     = xlim[2],
                  y0     = ylim[2],
                  x1     = xlim[2],
                  y1     = x$Summary[2, 1],
                  length = 0.1,
                  angle  = 30,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[2])
        } else {
            arrows(x0     = xlim[2],
                  y0     = x$Summary[2, 3],
                  x1     = xlim[2],
                  y1     = x$Summary[2, 1],
                  length = 0.1,
                  angle  = 90,
                  lwd    = lwd,
                  lty    = 1,
                  code   = 1,
                  col    = col[2])
        }
    }
    if(legend == TRUE) {
            legend(x      = xlim[1],
          y      = ylim[2],
          bty    = "n",
          legend = c("Genetic estimator",
                     "Field estimator"),
          pch    = pch,
          cex    = cex,
          lty    = lty,
          lwd    = lwd,
          col    = col)
    }
    # Reset margin settings
    par(mar = original_mar$mar)

}
