# updated 2017.03.19
# - added plot of reW and reL for ATE, using the code in making plots for the poster
#
# updated 2017.05.01
# - modified plot title, change color
#
# updated 2019.07.14
# - modified 



library(ggplot2)

# simulate rew and rel for single arm mean ----------------------------------------------------

# Figures here are included in the appendix of the paper.

rew <- function(q, py) {
    (0.5 * q * py + 1 - q)^{-1}
}

rel <- function(q, py_over_pl) {
    (q * py_over_pl + 1 - q)^{-1}
}



tmp <- data.frame(py = seq(0, 1, 0.01), q = 0.25)
tmp$re <- rew(tmp$q, tmp$py)
design_w <- tmp

tmp <- data.frame(py = seq(0, 1, 0.01), q = 0.1)
tmp$re <- rew(tmp$q, tmp$py)
design_w <- rbind(design_w, tmp)

design_w$q <- as.character(design_w$q)


tmp <- data.frame(py_over_pl = seq(0, 1, 0.01), q = 0.25)
tmp$re <- rel(tmp$q, tmp$py_over_pl)
design_l <- tmp

tmp <- data.frame(py_over_pl = seq(0, 1, 0.01), q = 0.1)
tmp$re <- rel(tmp$q, tmp$py_over_pl)
design_l <- rbind(design_l, tmp)

design_l$q <- as.character(design_l$q)


myggfont <- function(legend_text_size = 18,
                     legend_title_size = 20,
                     axis_text_size = 16,
                     axis_title_size = 20,
                     plot_title_size = 14) {
    ff <- theme(legend.text = element_text(size = legend_text_size),
                legend.title = element_text(size = legend_title_size),
                axis.text = element_text(size = axis_text_size),
                axis.title = element_text(size = axis_title_size, face="bold"),
                plot.title = element_text(size = plot_title_size))
    return(ff)
}

ggplot(design_w, aes(x = py, y = re, linetype = q)) + 
    geom_line() + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    # ggtitle("Estimating Mean Primary Outcome in One Arm") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_manual(guide = FALSE, values = c(2,1)) + 
    theme_bw() +
    myggfont() +
    # annotate("text", x = 0.25, y = c(1.32, 1.13), label = c(paste("{R['W;a']}^2 == 0.25"), paste("{R['W;a']}^2 == '0.10'")), size = 5, parse = TRUE) +
    # annotate("text", x = 0.40, y = c(1.32, 1.13), label = c(",", ","), size = 5) +
    # annotate("text", x = 0.53, y = c(1.32, 1.13), label = c(paste("{R['L|W;a']}^2 == 0"), paste("{R['L|W;a']}^2 == '0'")), size = 5, parse = TRUE)
    annotate("text", x = 0.26, y = c(1.32, 1.13), label = c(paste("{R['W;a']}^2 == 0.25"), paste("{R['W;a']}^2 == '0.10'")), size = 6, parse = TRUE) +
    annotate("text", x = 0.44, y = c(1.32, 1.13), label = c(",", ","), size = 5) +
    annotate("text", x = 0.61, y = c(1.32, 1.13), label = c(paste("{R['L|W;a']}^2 == 0"), paste("{R['L|W;a']}^2 == '0'")), size = 6, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW.png", width = 5, height = 5)


ggplot(design_l, aes(x = py_over_pl, y = re, linetype = q)) + 
    geom_line() + 
    xlab(expression(paste('p'[y],' / ', 'p'[l]))) +
    ylab("ARE") +
    # ggtitle("Estimating Mean Primary Outcome in One Arm") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_manual(guide = FALSE, values = c(2,1)) + 
    theme_bw() +
    myggfont() +
    # annotate("text", x = c(0.25, 0.25-0.15), y = c(1.32, 1.02), label = c(paste("{R['W;a']}^2 == 0"), paste("{R['W;a']}^2 == 0")), size = 5, parse = TRUE) +
    # annotate("text", x = c(0.36, 0.36-0.15), y = c(1.32, 1.02), label = c(",", ","), size = 5) +
    # annotate("text", x = c(0.53, 0.53-0.15), y = c(1.32, 1.02), label = c(paste("{R['L|W;a']}^2 == 0.25"), paste("{R['L|W;a']}^2 == '0.10'")), size = 5, parse = TRUE)
    annotate("text", x = c(0.24, 0.24-0.15), y = c(1.32, 1.02), label = c(paste("{R['W;a']}^2 == 0"), paste("{R['W;a']}^2 == 0")), size = 6, parse = TRUE) +
    annotate("text", x = c(0.37, 0.37-0.15), y = c(1.32, 1.02), label = c(",", ","), size = 5) +
    annotate("text", x = c(0.58, 0.58-0.15), y = c(1.32, 1.02), label = c(paste("{R['L|W;a']}^2 == 0.25"), paste("{R['L|W;a']}^2 == '0.10'")), size = 6, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REL.png", width = 5, height = 5)




# simulate rew and rel for estimating average treatment effect ----------------------------------------------------

# The following figures are included in Section 3 of the final paper.

rew <- function(q, py, gamma) {
    (0.5 * gamma * py + 1 - q)^{-1}
}

rel <- function(q, py_over_pl) {
    (q * py_over_pl + 1 - q)^{-1}
}


tmp <- expand.grid(py = seq(0, 1, 0.01), q = c(0.1, 0.25), gamma_multiplier = c(0, 0.5, 1))
tmp$gamma <- 2 * tmp$q * tmp$gamma_multiplier
tmp$re <- rew(tmp$q, tmp$py, tmp$gamma)
design_w <- tmp
design_w$q <- as.character(design_w$q)
design_w$gamma_multiplier <- as.character(design_w$gamma_multiplier)

tmp <- expand.grid(py_over_pl = seq(0, 1, 0.01), q = c(0.1, 0.25))
tmp$re <- rel(tmp$q, tmp$py_over_pl)
design_l <- tmp
design_l$q <- as.character(design_l$q)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#F0E442", "#0072B2", "#CC79A7")

ggplot(design_w, aes(x = py, y = re, linetype = q, color = gamma_multiplier)) + 
    geom_line() + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    # ggtitle("Estimating Average Treatment Effect") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_discrete(guide = FALSE) + 
    scale_color_manual(name = "Treatment effect\n heterogeneity", labels = c(
        bquote(paste(gamma, " = 0 (none)")),
        bquote(paste(gamma, " = ", {'R'[W]}^2, " (moderate)")),
        bquote(paste(gamma, " = 2", {'R'[W]}^2, " (maximal)"))),
        values = cbbPalette[1:3]) +
    theme_bw() +
    myggfont(plot_title_size = 17) +
    annotate("text", x = 0.15, y = c(1.36, 1.14), label = c(paste("{R[W]}^2 == 0.25"), paste("{R[W]}^2 == '0.10'")), size = 5, parse = TRUE) +
    annotate("text", x = 0.285, y = c(1.36, 1.14), label = c(",", ","), size = 5) +
    annotate("text", x = 0.40, y = c(1.36, 1.14), label = c(paste("{R['L|W']}^2 == 0"), paste("{R['L|W']}^2 == '0'")), size = 5, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW_ATE.png", width = 7.5, height = 5)

# revise the above figure so that it is all black and white, and only differentiated by line type.

design_w$mylinetype <- paste0(design_w$q, "_", design_w$gamma_multiplier)
ggplot(design_w, aes(x = py, y = re, linetype = mylinetype)) + 
    geom_line() + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    # ggtitle("Estimating Average Treatment Effect") +
    coord_cartesian(ylim = c(1, 1.4)) +
    # scale_linetype_discrete() +
    scale_linetype_manual(name = NULL, labels = c(
        bquote(paste({'R'[W]}^2, " = 0.1, ", gamma, " = 0")),
        bquote(paste({'R'[W]}^2, " = 0.1, ", gamma, " = ", {'R'[W]}^2)),
        bquote(paste({'R'[W]}^2, " = 0.1, ", gamma, " = 2", {'R'[W]}^2)),
        bquote(paste({'R'[W]}^2, " = 0.25, ", gamma, " = 0")),
        bquote(paste({'R'[W]}^2, " = 0.25, ", gamma, " = ", {'R'[W]}^2)),
        bquote(paste({'R'[W]}^2, " = 0.25, ", gamma, " = 2", {'R'[W]}^2))),
        values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")) +
        # values = 1:6) +
    theme_bw() +
    myggfont(plot_title_size = 17)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW_ATE_onlylinetype.png", width = 7.5, height = 5)


ggplot(design_w, aes(x = py, y = re, linetype = q, color = gamma_multiplier)) + 
    geom_line() + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_discrete(guide = FALSE) + 
    scale_color_manual(name = "Treatment effect\n heterogeneity", labels = c(
        bquote(paste(gamma, " = 0 : none")),
        bquote(paste(gamma, " = ", {'R'[W]}^2, " : moderate")),
        bquote(paste(gamma, " = 2", {'R'[W]}^2, " : maximal"))),
        values = cbbPalette[1:3]) +
    theme_bw() +
    myggfont() +
    annotate("text", x = 0.15, y = c(1.36, 1.14), label = c(paste("{R[W]}^2 == 0.25"), paste("{R[W]}^2 == '0.10'")), size = 5, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW_ATE_JSM2017.png", width = 7.5, height = 5)


ggplot(design_l, aes(x = py_over_pl, y = re, linetype = q)) + 
    geom_line() + 
    xlab(expression(paste('p'[y],' / ', 'p'[l]))) +
    ylab("ARE") +
    # ggtitle("Estimating Average Treatment Effect") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_manual(guide = FALSE, values = c(2,1)) + 
    theme_bw() +
    myggfont(plot_title_size = 17) +
    annotate("text", x = c(0.25, 0.25-0.15), y = c(1.32, 1.02), label = c(paste("{R['W']}^2 == 0"), paste("{R['W']}^2 == 0")), size = 5, parse = TRUE) +
    annotate("text", x = c(0.36, 0.36-0.15), y = c(1.32, 1.02), label = c(",", ","), size = 5) +
    annotate("text", x = c(0.53, 0.53-0.15), y = c(1.32, 1.02), label = c(paste("{R['L|W']}^2 == 0.25"), paste("{R['L|W']}^2 == '0.10'")), size = 5, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REL_ATE.png", width = 5, height = 5)


ggplot(design_l, aes(x = py_over_pl, y = re, linetype = q)) + 
    geom_line() + 
    xlab(expression(paste('p'[y],' / ', 'p'[l]))) +
    ylab("ARE") +
    coord_cartesian(ylim = c(1, 1.4)) +
    scale_linetype_manual(guide = FALSE, values = c(1,1)) + 
    theme_bw() +
    myggfont() +
    annotate("text", x = 0.25, y = c(1.32, 1.13), label = c(paste("{R[L]}^2 == 0.25"), paste("{R[L]}^2 == '0.10'")), size = 5, parse = TRUE)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REL_ATE_JSM2017.png", width = 5, height = 5)





# simulate rew and rel for JSM 2017 talk----------------------------------------------------



rew <- function(q, py, gamma) {
    (0.5 * gamma * py + 1 - q)^{-1}
}

rel <- function(q, py_over_pl) {
    (q * py_over_pl + 1 - q)^{-1}
}


tmp <- expand.grid(py = seq(0, 1, 0.01), q = 0.28)
tmp$gamma <- 0.06
tmp$re <- rew(tmp$q, tmp$py, tmp$gamma)
design_w <- tmp
design_w$q <- as.character(design_w$q)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#F0E442", "#0072B2", "#CC79A7")

ggplot(design_w, aes(x = py, y = re)) + 
    geom_line(size = 1.5) + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    ggtitle(expression(R[W]^2 == 0.28)) +
    coord_cartesian(ylim = c(1, 1.5)) +
    theme_bw() +
    myggfont() +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = 0.4, y = 1.31, label = paste("gamma == 0.06"), size = 7, parse = TRUE, color = "blue")
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW_ATE_1_JSM2017.png", width = 4, height = 4.3)


tmp <- expand.grid(py = seq(0, 1, 0.01), q = 0.28, gamma = c(0, 0.06, 2 * 0.28))
tmp$re <- rew(tmp$q, tmp$py, tmp$gamma)
design_w <- tmp
design_w$q <- as.character(design_w$q)
design_w$gamma <- factor(as.character(design_w$gamma))
design_w$gamma <- factor(design_w$gamma, levels(design_w$gamma)[c(2, 1, 3)])

ggplot(design_w, aes(x = py, y = re, linetype = gamma)) + 
    geom_line(size = 1.5) + 
    xlab(expression('p'[y])) +
    ylab("ARE") +
    ggtitle(expression(R[W]^2 == 0.28)) +
    scale_linetype_discrete(guide = FALSE) + 
    coord_cartesian(ylim = c(1, 1.5)) +
    theme_bw() +
    myggfont() +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = 0.4, y = c(1.1, 1.31, 1.42),
             label = c(paste("gamma == 0.56"), paste("gamma == 0.06"), paste("gamma == 0")), size = 7, parse = TRUE, color = "blue")
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REW_ATE_2_JSM2017.png", width = 4, height = 4.3)


tmp <- expand.grid(py_over_pl = seq(0, 1, 0.01), q = 0.28)
tmp$re <- rel(tmp$q, tmp$py_over_pl)
design_l <- tmp
design_l$q <- as.character(design_l$q)


ggplot(design_l, aes(x = py_over_pl, y = re)) + 
    geom_line(size = 1.5) + 
    xlab(expression(paste('p'[y],' / ', 'p'[l]))) +
    ylab("ARE") +
    ggtitle(expression(R["L|W"]^2 == 0.28)) +
    coord_cartesian(ylim = c(1, 1.5)) +
    theme_bw() +
    myggfont() +
    theme(plot.title = element_text(hjust = 0.5)) +
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/REL_ATE_JSM2017.png", width = 4, height = 4.3)



# simulate RedSS (reduction in sample size) ratio: r (for single arm mean) --------------------------------------------------------------

redss_ratio <- function(py, pl) {
    (1 - 0.5 * py) / (1 - pl^(-1) * py)
}

design <- expand.grid(py = seq(0, 1, 0.001), pl = seq(0, 1, 0.001))
design <- subset(design, py <= pl)

for (i in 1:nrow(design)) {
    if (i %% 1000 == 0) {
        print(i)
    }
    design$redss_ratio[i] <- redss_ratio(design$py[i], design$pl[i])
}

saveRDS(design, file = "~/Dropbox/Research/GST_TMLE/pl-py-redss_ratio.rds")
# saveRDS(design, file = "~/Dropbox/Research/GST_TMLE/pl-py-redss_ratio_smallset.rds")



# make plots 

design <- readRDS("~/Dropbox/Research/GST_TMLE/pl-py-redss_ratio.rds")

design$redss_ratio[design$redss_ratio >= 3] <- 3
design$redss_ratio[which(is.na(design$redss_ratio))] <- 1
names(design)[3] <- "r"




# make plot

myggfont <- function(legend_text_size = 18,
                     legend_title_size = 20,
                     axis_text_size = 16,
                     axis_title_size = 20,
                     plot_title_size = 20) {
    ff <- theme(legend.text = element_text(size = legend_text_size),
                legend.title = element_text(size = legend_title_size),
                axis.text = element_text(size = axis_text_size),
                axis.title = element_text(size = axis_title_size, face="bold"),
                plot.title = element_text(size = plot_title_size))
    return(ff)
}


ggplot(design) + 
    aes(x = pl, y = py, z = r, fill = r) + 
    geom_tile() + 
    xlab(expression('p'[l])) +
    ylab(expression('p'[y])) +
    ggtitle("Ratio of Sample Size Reductions\nfrom Adjusting for W and from Adjusting for L") +
    geom_contour(color = "white", breaks=c(1.2, 1.5, 2, 3)) + 
    coord_cartesian(xlim = c(0, 1.1)) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_fill_continuous(labels = c("  1.0", "  1.5", "  2.0", "  2.5", ">3.0")) + 
    theme_bw() +
    myggfont() +
    annotate("text", x = 1.08, y = c(0,0.3,0.5,0.67,0.8,1), label = c("r = 1", "r = 1.2", "r = 1.5", "r = 2.0", "r = 3.0", "r = \U221E"), size = 6)
ggsave("~/Dropbox/Research/GST_TMLE/writeup/figure/RedSS_ratio.png", width = 9, height = 8)
