## Evaluate model predict tissue PSI, per exon and per tissue
source("config.R")
tissues <- tissues_short

#########################################################
# Input:
#  - Predicted ASCOT test set with TSplice, run: 
#  - python evaluate.py -o data/test_ascot_performance.h5
#########################################################

# N <- 0 # at least vary N+1 tissue

library("cowplot")
library("ggrepel")

test_dt <- fread(paste0(DATA, "ascot_gtex_test.csv.gz"))
dim(test_dt)

measured <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "measured"))
pred <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "pred"))
pred_lm <- t(h5read(paste0(DATA, "ridge_ascot_test.h5"), "pred"))
x <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "x"))
dim(pred)
nna <- apply(measured, 2, function(x) sum(!is.na(x)))
rowm <- rowMeans(measured, na.rm = T)
vary <- abs(measured - rowm) > 0.2
nvary <- rowSums(vary, na.rm = T)

#
pcor <- sapply(c(1:dim(pred)[1]), function(i) cor(pred[i,], measured[i,], use="pairwise.complete.obs", m='p'))
scor <- sapply(c(1:dim(pred)[1]), function(i) cor(pred[i,], measured[i,], use="pairwise.complete.obs", m='s'))
scor_lm <- sapply(c(1:dim(pred_lm)[1]), function(i) cor(pred_lm[i,], measured[i,], use="pairwise.complete.obs", m='s'))
exon_rmse <- sapply(c(1:dim(pred)[1]), function(i) rmse(pred[i,], measured[i,]))
exon_cor = data.table(p=pcor, s=scor, rmse=exon_rmse)

# At least vary in 1 tissues
nna_exon <- apply(measured, 1, function(x) sum(!is.na(x)))
fter <- nvary > 0 & nna_exon >= 10
exon_cor <- exon_cor[fter]
exon_cor[, sum(s>0)] / dim(exon_cor)[1]

makdt <- function(i, fter=NULL){
  if(is.null(fter)){
    data.table(pred = pred[i,],
               measured = measured[i,],
               tissue = tissues)
  }else{
    data.table(pred = pred[fter,][i,],
               measured = measured[fter,][i,],
               tissue = tissues)
  }
}

example <- makdt(75)
example[, tissue := ifelse(pred > 0.2 | measured > 0.2, tissues, NA)]
exon_example <- ggplot(example, aes(pred, measured)) +
  geom_point() +
  labs(x=expression("Predicted "*Psi),
       y=expression("Measured "*Psi),
       title='Exon 9 of ABI2 gene') +
  geom_text_repel(aes(label=tissue), size=3, alpha=0.5) +
  xlim(0, 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")
exon_example
example[, rmse(pred, measured)]
example[, cor(pred, measured, m='s')]

hist_S <- ggplot(exon_cor, aes(s)) +
  geom_histogram() +
  labs(x = expression("Spearman correlation of predicted "*Psi*" per exon"))

exon_cor[, sum(s>0) / .N] 

## Compare to a base line with mean

dlogitpsi_mean <- logit(colMeans(measured, na.rm=T))
dlogitpsi_mean <- colMeans(x, na.rm=T)
base_line_exon <- sigmoid(
  sapply(
    logit(rowMeans(measured, na.rm = T)), function(i) {
      i + dlogitpsi_mean
    }
    )
  )
base_line_exon <- t(base_line_exon)

pcor_base <- sapply(c(1:dim(base_line_exon)[1]), 
               function(i) cor(base_line_exon[i,], measured[i,], use="pairwise.complete.obs", m='p'))
scor_base <- sapply(c(1:dim(base_line_exon)[1]), 
                    function(i) cor(base_line_exon[i,], measured[i,], use="pairwise.complete.obs", m='s'))
exon_rmse_base <- sapply(c(1:dim(base_line_exon)[1]), 
                         function(i) rmse(base_line_exon[i,], measured[i,]))
exon_cor_base = data.table(p_base=pcor_base, s_base=scor_base, s_lm=scor_lm, rmse_base=exon_rmse_base)
exon_cor_base <- exon_cor_base[fter]

exon_cor <- cbind(exon_cor, exon_cor_base)
exon_cor <- exon_cor[!is.na(s_base)]

ggplot(exon_cor, aes(s_lm, s)) +
  geom_point(size=0.5) +
  geom_abline(intercept = 0, slope = 1) +
  annotate('text', x=-0.6, y=0.7, 
           label=paste(exon_cor[, round(sum(s > s_lm) / .N, 3)]*100, '%')) +
  labs(title=expression('Correlation of per-exon predicted vs measured '*Psi),
       x='Ridge regression',
       y='TSplice') + 
  coord_cartesian(xlim = c(-0.7,0.9), ylim = c(-0.7,0.9))
ggsave(paste0(DATA, "baseline.pdf"), width = 6, height = 5)

exon_cor[, wilcox.test(s_base, s, paired = T)]

## Per tissue comparison
pred_rmse <- sapply(1:56, function(i) rmse(pred[fter,i], measured[fter,i]))
# rmse with mean
m <- rowMeans(measured, na.rm = T)
pred_mean <- sapply(1:56, function(i) rmse(m[fter], measured[fter,i]))

rmse_comp <- data.table(tissue=tissues, rmse_gain = pred_mean-pred_rmse)
rmse_comp[, tissue := factor(tissue, levels = tissues_short)]

performance <- ggplot(rmse_comp, aes(tissue, rmse_gain)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=-90, hjust = 0, vjust=0.5, size = 12)) +
  labs(x='', y='RMSE decrease')

up_row <- plot_grid(exon_example, hist_S, labels = c('A', 'B'),
                    align = 'h', rel_widths = c(1, 1.2), label_size = 18)
fig5 <- plot_grid(up_row, performance, labels = c('', 'C'),
                  ncol = 1, rel_heights = c(1, 1.3), label_size = 18)
ggsave('./data/Figure5.pdf', fig5, height=9, width = 10)





