### Plot the evaluation on predicting delta Logit PSI across tissues
source("config.R")

tissues <- tissues_short

#################################################################
# Input:
#  - Predicted ASCOT test set with TSplice, run: 
#  - python evaluate.py -o data/test_ascot_performance.h5
#################################################################

##' ## Load data
# measured PSI across tissues
measured <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "measured"))
# predicted delta logit PSI across tissues
x_hat <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "x_hat"))
x_hat_lm <- t(h5read(paste0(DATA, "ridge_ascot_test.h5"), "x_hat"))
# measured delta logit PSI across tissues
x <- t(h5read(paste0(DATA, "test_ascot_performance.h5"), "x"))
# Original exon table from ascot (test set)
exons <- fread(paste0(DATA, "ascot_gtex_test.csv.gz"))

##' ## Filter exons, at least vary in 1 tissue
rowm <- rowMeans(measured, na.rm = T)
# find tissues vary more than 0.2 from the mean
vary <- abs(measured - rowm) > 0.2
nvary <- rowSums(vary, na.rm = T)
nna_exon <- apply(measured, 1, function(x) sum(!is.na(x)))
fter <- nvary > 0 & nna_exon >= 10

# 
s_cor <- c()
s_cor_lm <- c()
p_cor <- c()
for(i in 1:56){
  s_cor <- c(s_cor, cor(c(x_hat[fter, i]), c(x[fter, i]), use="pairwise.complete.obs", m='s'))
  s_cor_lm <- c(s_cor_lm, cor(c(x_hat_lm[fter, i]), c(x[fter, i]), use="pairwise.complete.obs", m='s'))
  p_cor <- c(p_cor, cor(c(x_hat[fter, i]), c(x[fter, i]), use="pairwise.complete.obs", m='p'))
}

## compare with linear model TSplice baseline
lm_baseline <- qplot(s_cor_lm, s_cor) +
  geom_abline(intercept = 0, slope=1) +
  labs(y='TSplice', x='Ridge regression',
       title=expression('Correlation per-tissue predicted vs measured')) +
  xlim(-0.04, 0.42) + ylim(-0.04,0.42)
ggsave(filename=paste0(DATA, "Sup_lm_baseline.pdf"), width = 6, height = 5)

cor.test(s_cor, s_cor_lm) # p-value = 0.004102, cor=0.3777416 

# plot in ggplot
x <- as.data.table(x)
colnames(x) <- tissues
x[, exon_id := exons[, exon_id]]
x <- x[(fter)]
x <- melt(x, id.vars='exon_id', value.name='Measured', variable.name='tissues')
x_hat <- as.data.table(x_hat)
colnames(x_hat) <- tissues
x_hat[, exon_id := exons[, exon_id]]
x_hat <- x_hat[(fter)]
x_hat <- melt(x_hat, id.vars='exon_id', value.name="Predicted", variable.name='tissues')

merged <- merge(x, x_hat, by=c('exon_id', 'tissues'))

rho_labels <- data.table(tissues=tissues, 
                         rho=paste0("rho=", round(s_cor, 2)))

evaluate_dlogit <- ggplot(merged, aes(Predicted, Measured)) +
  geom_point(size=0.5) +
  geom_abline(slope=1, intercept=0, linetype = "dashed") +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  facet_wrap(~tissues, ncol=8) +
  geom_text(data=rho_labels, aes(x=3, y=10, label=rho)) +
  theme(strip.text.x = element_text(size=8))

# supplement Figure S1
ggsave(plot=evaluate_dlogit, filename=paste0(DATA, "evaluate_dlogit.pdf"), width = 12, height = 7)

rho <- data.table(tissues=tissues, 
                  rho=s_cor)
rho[, Brain := ifelse(substr(tissues, 1, 5) == 'Brain', "Brain", "Non-brain")]

# median example
example_plotx <- ggplot(merged[tissues=="Retina - Eye"], 
                        aes(Predicted, Measured)) +
  geom_point(size=0.5) + 
  ggtitle("Retina - Eye") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  annotate('text', 2.5, -9, label=paste("rho =", round(rho[tissues=='Retina - Eye', rho], 2)))

brain_nonbrain <- ggplot(rho, aes(Brain, rho)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  coord_flip() +
  theme(axis.title.y=element_blank())

cor_plot <- ggplot(rho, aes(tissues, rho)) +
  geom_bar(stat='identity') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5, size=10))

AC <- plot_grid(example_plotx, brain_nonbrain, labels = c('A', 'C'), ncol=1, label_size = 18)
plot_grid(AC, cor_plot, ncol=2, labels = c('', 'B'), rel_widths = c(1, 2.5))
ggsave(filename=paste0(DATA, "Figure4.pdf"), width = 10, height = 5)

