\name{error.dots}
\alias{error.dots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Show a  dot.chart with error bars for different groups or variables}
\description{Yet one more of the graphical ways of showing data with error bars for different groups.
A dot.chart with error bars for different groups or variables is found using from \code{\link{describe}},  \code{\link{describeBy}}, or  \code{\link{statsBy}}.
}
\usage{
error.dots(x, var = NULL, se = NULL, group = NULL, sd = FALSE, head = 12, tail = 12,
 sort = TRUE, decreasing = TRUE, main = NULL, alpha = 0.05, eyes = FALSE, min.n = NULL,
 max.labels = 40, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"),
  pt.cex = cex, pch = 21, gpch = 21, bg = par("bg"), color = par("fg"),
   gcolor = par("fg"), lcolor = "gray", xlab = NULL, ylab = NULL, xlim = NULL, ...)
}

\arguments{
  \item{x}{A data frame or matrix of raw data, or the resulting object from \code{\link{describe}},  \code{\link{describeBy}}, or  \code{\link{statsBy}} 

}
  \item{var}{The variable to show}

  \item{se}{Source of a standard error}


  \item{group}{A grouping variable, if desired}
  \item{sd}{if FALSE, confidence intervals in terms of standard errors, otherwise draw one standard deviation}
  \item{head}{The number of largest values to report}
  \item{tail}{The number of smallest values to report}     
  \item{sort}{Sort the groups/variables by value}
  \item{decreasing}{Should they be sorted in increasing or decreasing order (from top to bottom)}
  \item{main}{The caption for the figure}
  \item{alpha}{p value for confidence intervals}
  \item{eyes}{Draw catseyes for error limits}
  \item{min.n}{If using describeBy or statsBy, what should be the minimum sample size to draw}
  \item{max.labels}{Length of labels  (truncate after this value)}
  \item{labels}{Specify the labels versus find them from the row names}
  \item{groups}{ignored}

  \item{gdata}{ignored}

  \item{cex}{The standard meaning of cex for graphics}
  \item{pt.cex}{ignored}
  \item{pch}{Plot character}
  \item{gpch}{ignored}
  \item{bg}{background color}
  \item{color}{Color}
  \item{gcolor}{ignored}
  \item{lcolor}{ignored?}
  \item{xlab}{Label the x axis, if NULL, the variable name is used}
  \item{ylab}{If NULL, then the group rownames are used}
  \item{xlim}{If NULL, then calculated to show nice values}
  \item{\dots}{And any other graphic parameters we have forgotten}
}
\details{
Adapted from the dot.chart function to include error bars and to use the output of\code{\link{describe}},  \code{\link{describeBy}}, and  \code{\link{statsBy}}   
To speed up multiple plots, the function can work from the output of a previous run.  Thus describeBy will be done and the results can be show for multiple variables}
\value{
Returns (invisibily) either a describeBy or describe object
}
\references{Used in particular for showing http:\\sapa-project.org output.}
\author{William Revelle}


\seealso{
 \code{\link{describe}},  \code{\link{describeBy}}, or  \code{\link{statsBy}} as well as  \code{\link{error.bars}},  \code{\link{error.bars.by}}, or  \code{\link{statsBy}} 
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, var = NULL, se = NULL, group = NULL, sd = FALSE, 
    head = 12, tail = 12, sort = TRUE, decreasing = FALSE, main = NULL, 
    alpha = 0.05, eyes = FALSE, min.n = NULL, max.labels = 40, 
    labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), 
    gcolor = par("fg"), lcolor = "gray", xlab = NULL, ylab = NULL, 
    xlim = NULL, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    if (length(class(x)) > 1) {
        if (class(x)[1] == "psych") {
            obj <- class(x)[2]
            switch(obj, statsBy = {
                if (is.null(min.n)) {
                  se <- x$sd[, var]/sqrt(x$n[, var])
                  x <- x$mean[, var]
                } else {
                  se <- x$sd[, var]
                  n.obs <- x$n[, var]
                  x <- x$mean[, var]
                  if (sd) {
                    se <- x$sd[, var]
                  } else {
                    se <- se/sqrt(n.obs)
                  }
                  x <- subset(x, n.obs > min.n)
                  se <- subset(se, n.obs > min.n)
                  n.obs <- subset(n.obs, n.obs > min.n)
                }
            }, describe = {
                if (sd) {
                  se <- x$sd
                } else {
                  se <- x$se
                }
                labels <- rownames(x)
                x <- x$mean
                names(x) <- labels
            }, describeBy = {
                des <- x
                if (is.null(xlab)) xlab <- var
                var <- which(rownames(des[[1]]) == var)
                x <- se <- rep(NA, length(des))
                for (grp in 1:length(x)) {
                  x[grp] <- des[[grp]][["mean"]][var]
                  if (sd) {
                    se[grp] <- des[[grp]][["sd"]][var]
                  } else {
                    se[grp] <- des[[grp]][["se"]][var]
                  }
                }
                names(x) <- names(des)
                if (is.null(xlab)) xlab <- var
            })
        }
    }
    else {
        if (is.null(group)) {
            des <- describe(x)
            x <- des$mean
            if (sd) {
                se <- des$sd
            }
            else {
                se <- des$se
            }
            names(x) <- rownames(des)
        }
        else {
            if (is.null(xlab)) 
                xlab <- var
            des <- describeBy(x, group = group)
            x <- se <- rep(NA, length(des))
            names(x) <- names(des)
            var <- which(rownames(des[[1]]) == var)
            for (grp in 1:length(des)) {
                x[grp] <- des[[grp]][["mean"]][var]
                if (sd) {
                  se[grp] <- des[[grp]][["sd"]][var]
                }
                else {
                  se[grp] <- des[[grp]][["se"]][var]
                }
            }
        }
    }
    n.var <- length(x)
    if (sort) {
        ord <- order(x, decreasing = decreasing)
    }
    else {
        ord <- n.var:1
    }
    x <- x[ord]
    se <- se[ord]
    temp <- temp.se <- rep(NA, min(head + tail, n.var))
    if ((head + tail) < n.var) {
        if (head > 0) {
            temp[1:head] <- x[1:head]
            temp.se[1:head] <- se[1:head]
            names(temp) <- names(x)[1:head]
        }
        if (tail > 0) {
            temp[(head + 1):(head + tail)] <- x[(length(x) - 
                tail + 1):length(x)]
            temp.se[(head + 1):(head + tail)] <- se[(length(x) - 
                tail + 1):length(x)]
            names(temp)[(head + 1):(head + tail)] <- names(x)[(length(x) - 
                tail + 1):length(x)]
        }
        x <- temp
        se <- temp.se
    }
    if (missing(main)) {
        if (sd) {
            main <- "means + standard deviation"
        }
        else {
            main = "Confidence Intervals around the mean"
        }
    }
    labels <- names(x)
    if (sd) {
        ci <- se
    }
    else {
        ci <- qnorm((1 - alpha/2)) * se
    }
    if (!is.null(ci) && is.null(xlim)) 
        xlim <- c(min(x - ci), max(x + ci))
    labels <- substr(labels, 1, max.labels)
    if (eyes) {
        ln <- seq(-3, 3, 0.1)
        rev <- (length(ln):1)
    }
    if (!is.numeric(x)) 
        stop("'x' must be a numeric vector or matrix")
    n <- length(x)
    if (is.matrix(x)) {
        if (is.null(labels)) 
            labels <- rownames(x)
        if (is.null(labels)) 
            labels <- as.character(1L:nrow(x))
        labels <- rep_len(labels, n)
        if (is.null(groups)) 
            groups <- col(x, as.factor = TRUE)
        glabels <- levels(groups)
    }
    else {
        if (is.null(labels)) 
            labels <- names(x)
        glabels <- if (!is.null(groups)) 
            levels(groups)
        if (!is.vector(x)) {
            warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
            x <- as.numeric(x)
        }
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1L:n
        y <- o
        ylim <- c(0, n + 1)
    }
    else {
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
            col = color, las = 2, cex = cex, ...)
    }
    abline(h = y, lty = "dotted", col = lcolor)
    points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)
    if (!is.null(ci)) {
        if (!eyes) {
            segments(x - ci, y, x + ci, y, col = par("fg"), lty = par("lty"), 
                lwd = par("lwd"))
        }
        else {
            catseyes(x, y, se = se, n = NULL, alpha = alpha, 
                density = -10)
        }
    }
    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                cex = pt.cex/cex, ...)
        }
    }
    axis(1)
    box()
    invisible()
    if (!is.null(group)) 
        result <- des
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
