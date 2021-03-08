source("config.R")
tissues <- tissues_short

#################################################
# Figure - Tissue specific splicing, description
# Input: 
#   - All ASCOT cassette exon PSI per tissue
#################################################
library("cowplot")
# all data
ref_psi <- h5read(paste0(DATA, "ref_psi_all.h5"), "ref_psi") %>% t
rowm <- rowMeans(ref_psi, na.rm = T)
# find tissues vary more than 0.1 from the mean
vary <- abs(ref_psi - rowm) > 0.1
nvary <- rowSums(vary, na.rm = T)
tb <- nvary %>% table
# barplot(tb[2:36])
tb <- as.data.table(tb)
colnames(tb) <- c('N', 'count')
tb[, N := as.numeric(N)]
tb[, g := N > 9]

hist_as <- ggplot(tb[2:57], aes(N, count, fill=g)) +
  geom_bar(stat='identity') +
  coord_cartesian(xlim = c(1,56)) + 
  scale_x_continuous(breaks = c(1,10,20,30, 40)) +
  labs(x='Differentially spliced tissues (+/-10% from average PSI)',
       y='Exons') +
  theme(legend.position = "none", plot.margin = unit(c(0.5,4,0.5,4), "cm")) +
  scale_fill_manual(values=c("#999999", "#56B4E9"))
hist_as
# fraction of differentially spliced exons
sum(nvary >= 1) / length(nvary)
sum(nvary > 9)

# heatmap: vary in at least 10 tissues
## First impute to do clustering
imputed_m <- ref_psi[nvary > 9,] # impute rowmeans
k <- which(is.na(imputed_m), arr.ind=TRUE)
imputed_m[k] <- rowMeans(imputed_m, na.rm=TRUE)[k[,1]]
# Then plot with the none imputed data
mdt <- data.table(ref_psi)
colnames(mdt) <- tissues
mdt <- mdt[nvary > 9] # at least vary in 10 tissues
mdt[, exons := seq(dim(mdt)[1])]
mdt <- melt(mdt, id.vars = 'exons', variable.name ='Tissues', value.name = 'PSI')

#data <- scale(t(mdt))
hc <- hclust( dist(imputed_m, method = "euclidean"), method = "ward.D" )
mdt[, exons := factor(exons, levels=hc$order)]
# library(ggdendro)
# ggdendrogram(hc, labels = FALSE, rotate = 90)

ordx <- hclust( dist(t(imputed_m), method = "euclidean"), method = "ward.D" )$order
mdt[, Tissues := factor(Tissues, levels=tissues[ordx])]

ref_map <- ggplot(mdt, aes(Tissues, exons)) +
  geom_tile(aes(fill = PSI)) +
  scale_fill_gradient2(low = 'blue',
                       mid='white',
                       high = 'red',
                       midpoint=0.5,
                       space ="Lab") +
  labs(y = 'Exons') +
  theme(axis.ticks = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle = -90, hjust=0, vjust=0.5, size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9),
        legend.key.size = unit(0.2,"cm"),
        legend.text=element_text(size=9),
        legend.title=element_text(size=9))

fig_heatmap <- plot_grid(hist_as, ref_map, labels=c('A','B'), ncol=1, rel_heights = c(0.7,2))
ggsave(paste0(DATA, "Figure1.pdf"), fig_heatmap, width = 7, height = 12)

#################################################
# Figure - Tissue specific variant effect, description

# Input generated in:
#  - GTEx-variant-effect-prediction.ipynb

#################################################
ref_psi <- h5read(paste0(DATA, "Variant_benchmark.h5"), "ref_psi") %>% t
alt_psi <- h5read(paste0(DATA, "Variant_benchmark.h5"), "alt_psi") %>% t
homo <- h5read(paste0(DATA, "Variant_benchmark.h5"), "homo")
homo <- homo == 1
ref_psi <- ref_psi[homo, ]
alt_psi <- alt_psi[homo, ]
dim(ref_psi)  #  1767   53

dpsi <- alt_psi - ref_psi
mdpsi <- rowMeans(dpsi, na.rm = T)
vary <- abs(dpsi - mdpsi) > 0.2
nvary <- rowSums(vary, na.rm = T)
sum(vary, na.rm = TRUE)  # 1030

dpsi_variation <- data.frame(mdpsi = rep(mdpsi, 53), dpsi = c(dpsi))
mdpsi_plot <- ggplot(data=dpsi_variation, aes(mdpsi, dpsi)) +
  labs(x=expression('Mean '*Delta*Psi* ' across tissues'),
       y=expression('Measured '*Delta*Psi*' (all tissues)')) +
  theme(legend.position = 'none') +
  geom_hex(bins = 70) + 
  #geom_point(alpha=0.5) +
  #scale_colour_brewer(palette = "Greys")
  scale_fill_distiller(name = "count", trans = "log10", palette = "GnBu")

ggsave(paste0(DATA, "Figure1_as_mdpsi_scatter.pdf"), width = 5, height = 3.5)

ecdf_dt <- data.table(x=abs(rep(mdpsi, 53) - c(dpsi)))
ecdf_dt <- ecdf_dt[!is.na(x)]
ecdf_dt[, len := .N]
ecdf_fn <- ecdf(ecdf_dt[,x])
highlight2 <- 1-ecdf_fn(0.2) %>% round(3)
ecdf_dpsi <- ggplot(ecdf_dt, aes(x)) +
  geom_step(stat='ecdf') + 
  scale_x_reverse() +
  scale_y_continuous(breaks = c(0, highlight2, 0.25, 0.5, 0.75, 1)) +
  background_grid() +
  geom_hline(yintercept = highlight2, linetype=2) +
  labs(x=expression('|'*Delta*Psi*'-mean('*Delta*Psi*')|'), y='quantile')


plt <- plot_grid(mdpsi_plot, ecdf_dpsi, labels = c('A','B'))
ggsave(filename = paste0(DATA, "Figure2.pdf"), width = 10, height = 5, plot = plt)
