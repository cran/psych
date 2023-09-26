"dia.arrow" <- function (from, to, labels = NULL, scale = 1, cex = 1, adj = 2, 
    both = FALSE, pos = NULL, l.cex, gap.size = NULL, draw = TRUE, 
    col = "black", lty = "solid", ...) 
{
    if (missing(gap.size)) 
        gap.size <- 0.2
    if (missing(l.cex)) 
        l.cex <- cex
    text.values <- NULL
    radius1 <- radius2 <- 0
    if (!is.list(to)) {
        tocenter <- to
    }
    else {
        tocenter <- to$center
    }
    if (is.list(from)) {
        if (!is.null(from$radius)) {
            radius1 <- from$radius
            radius2 <- 0
        }
        if (!is.null(from$center)) 
            from <- from$center
    }
    if (is.list(to)) {
        if (!is.null(to$radius)) {
            radius2 <- to$radius
            to <- to$center
        }
    }
    theta <- atan((from[2] - to[2])/(from[1] - to[1]))
    costheta <- cos(theta)
    sintheta <- sin(theta)
    dist <- sqrt((to[1] - from[1])^2 + (to[2] - from[2])^2)
    if ((adj > 3) || (adj < 1)) {
        x <- (to[1] + from[1])/2
        y <- (to[2] + from[2])/2
    }
    else {
        x <- from[1] - sign(from[1] - to[1] + 0.001) * (4 - adj) * 
            costheta * dist/4
        y <- from[2] - sign(from[1] - to[1] + 0.001) * (4 - adj) * 
            sintheta * dist/4
    }
    if (is.null(labels)) {
        h.size <- 0
    }
    else {
        h.size <- nchar(labels) * cex * gap.size
    }
    if (is.null(labels)) {
        v.size <- 0
    }
    else {
        v.size <- cex * 0.7 * gap.size
    }
    if (from[1] < to[1]) {
        h.size <- -h.size
        radius1 <- -radius1
        radius2 <- -radius2
    }
    x0 <- from[1] - costheta * radius1
    xr <- x + h.size * costheta * v.size
    y0 <- from[2] - sintheta * radius1
    xl <- x - h.size * costheta * v.size
    yr <- y + v.size * sintheta * h.size
    yl <- y - v.size * sintheta * h.size
    xe <- to[1] + costheta * radius2
    ye <- to[2] + sintheta * radius2
    if (!is.null(labels)) {
        if (draw) {
            text(x, y, labels, cex = l.cex, pos = pos, ...)
            text.values <- NULL
        }
        else {
            if (is.null(pos)) 
                pos <- 0
            text.values <- c(x, y, labels, pos = pos, cex = l.cex)
        }
    }
    if (draw) {
        arrows(x0, y0, xr, yr, length = (both + 0) * 0.1 * scale, 
            angle = 30, code = 1, col = col, lty = lty, ...)
        arrows(xl, yl, xe, ye, length = 0.1 * scale, angle = 30, 
            code = 2, col = col, lty = lty, ...)
        arrow.values <- NULL
    }
    else {
        arrow.values <- list(x0, y0, xr, yr, length = (both + 
            0) * 0.1 * scale, angle = 30, code = 1, xl, yl, xe, 
            ye, length2 = 0.1 * scale, angle2 = 30, code2 = 2, 
            col = col, lty = lty, ...)
    }
    invisible(list(text.values = text.values, arrow.values = arrow.values))
}
