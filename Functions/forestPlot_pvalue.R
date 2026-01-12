forestPlot_pvalue <- function (mafCompareRes,
                        pVal = 0.05,
                        fdr = NULL,
                        color = c("maroon","royalblue"),
                        geneFontSize = 0.8,
                        titleSize = 1.2,
                        lineWidth = 1,
                        pval_digits = 2,
                        pval_scientific = TRUE
)
{
  if (is.null(fdr)) {
    m.sigs <- mafCompareRes$results[pval < pVal]
  } else {
    m.sigs <- mafCompareRes$results[adjPval < fdr]
  }
  
  m1Name <- mafCompareRes$SampleSummary[1, Cohort]
  m2Name <- mafCompareRes$SampleSummary[2, Cohort]
  m1.sampleSize <- mafCompareRes$SampleSummary[1, SampleSize]
  m2.sampleSize <- mafCompareRes$SampleSummary[2, SampleSize]
  
  if (length(color) == 1) {
    vc_col <- c(color, color)
  } else if (length(color) == 2) {
    vc_col <- color
  } else {
    stop("colors length must be less equal than 2")
  }
  
  if (is.null(names(vc_col))) {
    names(vc_col) <- c(m2Name, m1Name)
  } else if (!(m1Name %in% names(vc_col)) && !(m2Name %in% names(vc_col))) {
    stop(paste0("\ninput named vector [color] must contain both group name, \nwhich should be like c(",
                m1Name, ", ", m2Name, ") but now are ", list(names(vc_col)), "\n"))
  } else if (!(m1Name %in% names(vc_col)) || !(m2Name %in% names(vc_col))) {
    stop(paste0("\nif you pass [color] with named vector, then both of the names matching group names must be Explicitly declared,\n",
                "which should be like c(", m1Name, ", ", m2Name, ") but now are ", list(names(vc_col)), "\n"))
  }
  
  if (nrow(m.sigs) < 1) stop("No differetially mutated genes found !")
  
  m.sigs$Hugo_Symbol <- factor(x = m.sigs$Hugo_Symbol, levels = rev(m.sigs$Hugo_Symbol))
  m.sigs$or_new <- ifelse(m.sigs$or > 3, 3, m.sigs$or)
  m.sigs$upper <- ifelse(m.sigs$ci.up  > 3, 3, m.sigs$ci.up)
  m.sigs$lower <- ifelse(m.sigs$ci.low > 3, 3, m.sigs$ci.low)
  m.sigs$pos <- rev(1:nrow(m.sigs))
  m.sigs <- m.sigs[order(pos)]
  m.sigs$pos <- 1:nrow(m.sigs)

  .format_p_smart <- function(x, dp = 3, sig = 3) {
    x <- as.numeric(x)
    out <- formatC(x, format = "f", digits = dp)
    is_zero <- suppressWarnings(as.numeric(out) == 0 & !is.na(x))
    if (any(is_zero)) {
      out[is_zero] <- format(signif(x[is_zero], sig), scientific = FALSE, trim = TRUE)
    }
    out
  }
  
  m.sigs$p_label <- .format_p_smart(m.sigs$pval, dp = 3, sig = 3)
  if ("adjPval" %in% names(m.sigs)) {
    m.sigs$padj_label <- .format_p_smart(m.sigs$adjPval, dp = 3, sig = 3)
  }
  
  xlims <- c(0, 4)
  ylims <- c(0.75, nrow(m.sigs))
  graphics::layout(
    mat = matrix(c(1,2,3,4,5, 6,6,6,6,6), byrow = TRUE, ncol = 5, nrow = 2),
    widths  = c(3.8, 1.0, 1.0, 1.0, 1.6),
    heights = c(6, 1.2)
  )
  
  par(mar = c(3, 0.4, 3, 2), xaxs = "i")
  plot(NA, xlim = c(0, 4), ylim = ylims, axes = FALSE, xlab = NA, ylab = NA)
  
  apply(m.sigs[, .(or, ci.up, ci.low, ci.up, or_new, upper, lower, pos)], 1, function(x) {
    p <- x[5]
    u <- x[6]; u_orig <- x[2]
    l <- x[7]; l_orig <- x[3]
    ypos <- x[8]
    
    if (p < 1) {
      linecolor <- vc_col[m2Name]
    } else if (p > 1) {
      linecolor <- vc_col[m1Name]
    } else {
      linecolor <- "black"
    }
    
    points(x = p, y = ypos, pch = 16, cex = 1.1 * (lineWidth))
    segments(x0 = l, y0 = ypos, x1 = u, y1 = ypos, lwd = lineWidth, col = linecolor)
    
    if (u_orig > 3) {
      segments(x0 = 3, y0 = ypos, x1 = 3.25, y1 = ypos, lwd = lineWidth, col = linecolor)
      points(x = 3.25, y = ypos, pch = ">", cex = 1.1 * (lineWidth))
    }
    if (l_orig > 3) {
      segments(x0 = 3, y0 = ypos, x1 = 3.25, y1 = ypos, lwd = lineWidth, col = linecolor)
      points(x = 3.25, y = ypos, pch = ">", cex = 1.1 * (lineWidth))
    }
  })
  
  abline(v = 1, lty = 2, col = "gray", xpd = FALSE)
  axis(side = 1, at = 0:3, labels = c(0:3), font = 1, pos = 0.5, cex.axis = 1.3)
  mtext(text = m.sigs$Hugo_Symbol, side = 4, line = 0.2, at = 1:nrow(m.sigs),
        font = 3, las = 2, cex = geneFontSize, adj = 0)
  
  mtitle <- paste(m2Name, " (n = ", m2.sampleSize, ")", " v/s ",
                  m1Name, " (n = ", m1.sampleSize, ")", sep = "")
  title(main = mtitle, font = 1, adj = 0, cex.main = titleSize)
  
  par(mar = c(3, 0, 3, 0))
  plot(NA, xlim = c(0, 1), ylim = ylims, axes = FALSE, xlab = "", ylab = "")
  text(x = 0.5, y = 1:nrow(m.sigs),
       labels = as.numeric(unlist(m.sigs[, 3])),
       adj = 0.5, cex = 1.2 * geneFontSize)
  title(main = m2Name, cex.main = min(1, titleSize), line = 1)
  
  par(mar = c(3, 0, 3, 0))
  plot(NA, xlim = c(0, 1), ylim = ylims, axes = FALSE, xlab = "", ylab = "")
  text(x = 0.5, y = 1:nrow(m.sigs),
       labels = as.numeric(unlist(m.sigs[, 2])),
       adj = 0.5, cex = 1.2 * geneFontSize)
  title(main = m1Name, cex.main = min(1, titleSize), line = 1)
  
  
  par(mar = c(3, 0, 3, 0))
  plot(NA, xlim = c(0, 1), ylim = ylims, axes = FALSE, xlab = "", ylab = "")
  text(x = 0.5, y = 1:nrow(m.sigs),
       labels = ifelse(is.finite(m.sigs$or), sprintf("%.3f", m.sigs$or), "Inf"),
       adj = 0.5, cex = 1.2 * geneFontSize)
  title(main = "OR", cex.main = min(1, titleSize), line = 1)
  
  m.sigs$significance <- ifelse(as.numeric(m.sigs$pval) < 0.001, "***",
                                ifelse(as.numeric(m.sigs$pval) < 0.01,  "**",
                                       ifelse(as.numeric(m.sigs$pval) < 0.05,  "*", "NS")))
  
  par(mar = c(3, 0, 3, 0))
  plot(NA, xlim = c(0, 1), ylim = ylims, axes = FALSE, xlab = "", ylab = "")
  text(x = 0.5, y = 1:nrow(m.sigs),
       labels = paste0(m.sigs$significance, "  (p=", m.sigs$p_label, ")"),
       adj = 0.5, cex = 1.2 * geneFontSize) 
  title(main = "P-value", cex.main = min(1, titleSize), line = 1)
  
  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = NA, ylab = NA)
  text(x = 0, labels = paste0("Odds ratio with 95% CI\n(1 = no effect, < 1 ",
                              m2Name, " has more mutants)"), y = 0.6, adj = 0, xpd = TRUE, cex = 1.2)
}