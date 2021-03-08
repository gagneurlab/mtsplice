source("config.R")
library(cowplot)
library(ggrepel)
tissues <- tissues_gtex

#################################################
# Figure - GTEx variant effect prediction benchmark

# Input generated in:
#  - GTEx-variant-effect-prediction.ipynb

#################################################

benchmark <- "Variant_benchmark.h5"

ref_psi <- t(h5read(paste0(DATA, benchmark), "ref_psi"))
alt_psi <- t(h5read(paste0(DATA, benchmark), "alt_psi"))
MMSplice_MTSplice <- t(h5read(paste0(DATA, benchmark), "MMSplice_MTSplice"))
MMSplice <- t(h5read(paste0(DATA, benchmark), "MMSplice"))
dPSI_Measured <- t(h5read(paste0(DATA, benchmark), "dPSI_Measured"))

## 
vary <- dPSI_Measured - rowMeans(dPSI_Measured, na.rm = T)
nvary <- apply(vary, 1, function(x) sum(abs(x) > 0.2, na.rm = T))
selected <- nvary > 0
sum(selected)
vary_instances <- abs(vary) > 0.2
vary_instances[is.na(vary_instances)] <- FALSE
vary_pred_mm <- MMSplice - rowMeans(MMSplice, na.rm = T)

#' Evaluate per tissue
# plot dPSI scatter plot per tissue for MMSplice and MTSplice

dpsi_plot_tissue <- function(i, select=selected, title='Amygdala'){
  dt <- data.table(mmsplice = MMSplice[select,i],
                   mtsplice = MMSplice_MTSplice[select,i],
                   dpsi = dPSI_Measured[select,i],
                   vary = selected[select]) %>% 
    melt(id.vars=c('dpsi', 'vary'))
  dt[, cor(dpsi, value, use='pairwise.complete.obs'), by=variable]
  dt[, rmse(dpsi, value), by=variable]
  labels <- c(mmsplice = "MMSplice", mtsplice = "MTSplice")
  mselabel <- dt[, rmse(dpsi, value), by=variable][, V1]
  mselabel <- sapply(mselabel, format_mse)
  label_dt <- data.frame(value=-0.5, dpsi=0.3, label=mselabel[1])
  mm_dt <- dt[variable=='mmsplice']
  plt_mm <- ggplot(mm_dt, aes(value, dpsi)) +
    geom_point(aes(color=vary)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed') +
    geom_text(aes(label=label), data=label_dt) +
    labs(x = expression("Predicted "*Delta*Psi), 
         y = expression("Measured "*Delta*Psi),
         title = paste('MMSplice -', title)) +
    theme(strip.background = element_rect(colour="white", fill="white")) +
    theme(legend.position="bottom")
  
  mm_dt <- dt[variable=='mtsplice']
  label_dt <- data.frame(value=-0.5, dpsi=0.3, label=mselabel[2])
  plt_mt <- ggplot(mm_dt, aes(value, dpsi)) +
    geom_point(aes(color=vary)) +
    geom_abline(intercept = 0, slope = 1, linetype='dashed') +
    geom_text(aes(label=label), data=label_dt) +
    labs(x = expression("Predicted "*Delta*Psi), 
         y = expression("Measured "*Delta*Psi),
         title = paste('MTSplice -', title)) +
    theme(strip.background = element_rect(colour="white", fill="white")) +
    theme(legend.position="bottom")
  plot_grid(plt_mm, plt_mt, labels=c('A','B'))
}

dpsi_tissue_compare <- dpsi_plot_tissue(10,  rep(TRUE, dim(MMSplice)[1]))
dpsi_tissue_compare

### Compare RMSE for selected variants
mse_compare <- sapply(seq(dim(MMSplice)[2]), function(i){
  c(rmse(MMSplice[selected,i], dPSI_Measured[selected,i]),
    rmse(MMSplice_MTSplice[selected,i], dPSI_Measured[selected,i]))
})
mse_compare <- data.table(t(mse_compare))
colnames(mse_compare) <- c("MMSplice", "MMSplice_MTSplice")
mse_compare[, Brain := c(rep('Non-Brain', 8), rep('Brain', 12), rep('Non-Brain', 33))]

mse_compare[10]  # example tissue

testis <- mse_compare[MMSplice - MMSplice_MTSplice < -0.01]
testis[, Val := c(-0.0065)]
testis[, label := c('triangle down filled')]

all_tissues <- ggplot(mse_compare, aes(Brain, MMSplice - MMSplice_MTSplice)) +
  geom_boxplot(outlier.color = 'white') +
  geom_jitter(width = 0.3, size=3) +
  labs(x='', y='RMSE(MMSplice) - RMSE(MTSplice)') +
  coord_cartesian(ylim = c(-0.006, 0.01)) +
  geom_hline(yintercept = 0, linetype='dashed') +
  annotate('text', x=0, y=0.05, label='S') +
  geom_point(data=testis, aes(Brain, Val, shape=label), size=4) +
  scale_shape_identity() + 
  geom_text(data=testis, aes(Brain, Val+0.0008, label='Testis \n -0.0374'), size=4)

# compare only tissue large vary
mse_compare_vary <- sapply(seq(dim(MMSplice)[2]), function(i){
  sl <- vary_instances[,i]
  c(rmse(MMSplice[sl,i], dPSI_Measured[sl,i]),
    rmse(MMSplice_MTSplice[sl,i], dPSI_Measured[sl,i]))
})

diff_vary <- mse_compare_vary[1,] - mse_compare_vary[2,]
mse_compare_vary <- data.table(t(mse_compare_vary))
colnames(mse_compare_vary) <- c("MMSplice", "MMSplice_MTSplice")
label_text <- c("Testis")
mse_compare_vary[, label := tissues_gtex]
mse_compare_vary[, label := ifelse(label %in% label_text, label, NA)]
mse_compare_vary[, Brain := c(rep('Non-Brain', 8), rep('Brain', 12), rep('Non-Brain', 33))]
mse_compare_vary <- mse_compare_vary[apply(vary_instances, 2, sum) > 10]

# 0.2444            0.2405                           Testis Non-Brain

vary_tissues <- ggplot(mse_compare_vary, aes(Brain, MMSplice - MMSplice_MTSplice)) +
  geom_boxplot(outlier.color = 'white') +
  geom_jitter(width = 0.3, size=3) +
  labs(x='', y='RMSE(MMSplice) - RMSE(MTSplice)') +
  geom_text(aes(label = label), na.rm = TRUE, hjust = -0.5) +
  geom_hline(yintercept = 0, linetype='dashed')

gtex_tissue_specific <- plot_grid(all_tissues, vary_tissues, labels = c('A','B'))
ggsave(paste0(DATA, "Sup_GTEx.pdf"), plot = gtex_tissue_specific, width = 12.5, height = 5.5)

# remove tissues with less than 10 data points
valid_tissues <- apply(dPSI_Measured, 2, function(x) sum(!is.na(x))) > 10

# label tissue MTSplice is worse
mse_compare[, text_label := ifelse(mse_compare[, abs(MMSplice_MTSplice - MMSplice) > 0.002], tissues_gtex, NA)]
mse_compare[, text_label := ifelse(text_label=='Brain - Spinal cord (cervical c-1)', "Brain - Spinal cord", text_label)]
mse_plot <- ggplot(mse_compare[valid_tissues], aes(MMSplice, MMSplice_MTSplice)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype='dashed') +
  labs(y = 'MTSplice', title='Root mean square error') +
  geom_text_repel(aes(label=text_label), size=4, alpha=0.5, 
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"),
                  nudge_x       = 0.002) +
  xlim(0.12, 0.19) + ylim(0.12, 0.2)

mse_compare[valid_tissues][, wilcox.test(MMSplice, MMSplice_MTSplice, paired = T)]
mse_compare[valid_tissues][, sum(MMSplice - MMSplice_MTSplice > 0)]

# All variants 1767
mse_compare_all <- sapply(seq(dim(MMSplice)[2]), function(i){
  c(rmse(MMSplice[,i], dPSI_Measured[,i]),
    rmse(MMSplice_MTSplice[,i], dPSI_Measured[,i]))
})
mse_compare_all <- data.table(t(mse_compare_all))
colnames(mse_compare_all) <- c("MMSplice", "MMSplice_MTSplice")

mse_compare_all[10]  # example tissue
mse_compare_all[valid_tissues][, wilcox.test(MMSplice, MMSplice_MTSplice, paired = T)]
mse_compare_all[valid_tissues][, sum(MMSplice - MMSplice_MTSplice > 0)]

tissue_evaluation <- plot_grid(dpsi_tissue_compare, mse_plot, nrow = 1, rel_widths = c(2,1), labels=c('', 'C'))
ggsave(paste0(DATA, "Figure6.pdf"), plot = tissue_evaluation, width = 14, height = 4.5)
