# Tissue-specific variant effect prediction revealed prioritize splicing variants of autism disease
source("config.R")
library(cowplot)
library(ggrepel)

#################################################
# Figure - Application of MTSplice to Autism

# Input generated in:
#  - notebooks/Autism/autism-MMSplice-MTSplice.ipynb
#  - notebooks/Autism/autism-MMSplice.ipynb
#################################################

# Tissue-specific variant effect prediction revealed prioritize splicing variants of autism disease
tissue_all_variants <- fread(paste0(DATA, "MMSplice_pred_oe.csv"), drop=1)
tissue_all_variants <- tissue_all_variants[abs(delta_logit_psi) > 0.05]
dim(tissue_all_variants)
tissue_all_variants[, sum(Proband==T)]
tissue_all_variants[, sum(Proband==F)]

all_p <- wilcox.test(tissue_all_variants[Proband==T, delta_logit_psi],
                     tissue_all_variants[Proband==F, delta_logit_psi],
                     alternative='less')$p.value

lof_p <- wilcox.test(tissue_all_variants[(lof_gene)][Proband==T, delta_logit_psi],
                     tissue_all_variants[(lof_gene)][Proband==F, delta_logit_psi], 
                     alternative='less')$p.value

p.adjust(c(all_p, lof_p), method='BH')

all_variants_mmsplice <- summarySE(tissue_all_variants, 
                                   measurevar = 'delta_logit_psi', 
                                   groupvars = c('Proband'))
lof_variants_mmsplice <- summarySE(tissue_all_variants[(lof_gene)], 
                                   measurevar = 'delta_logit_psi', 
                                   groupvars = c('Proband'))
variants_mmsplice <- rbindlist(list(all_variants_mmsplice, 
                                    lof_variants_mmsplice), 
                               idcol = 'group')
variants_mmsplice[, group := as.factor(group)]
variants_mmsplice[, Proband := ifelse(Proband, 'Proband', 'Sibling')]

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

mmsplice <- ggplot(variants_mmsplice, aes(group, delta_logit_psi, color=Proband)) +
  geom_pointrange(aes(ymin=delta_logit_psi-ci, ymax=delta_logit_psi+ci), 
                  position = position_dodge(width = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  labs(x=element_blank(), y=expression('MMSplice '*Delta*'logit'*Psi)) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = 'bottom') +
  #ylim(-0.3,-0.05) +
  annotate("text", x=1, y=-0.05, label=paste("italic(P) == ", format.pval(all_p, digits = 2)), parse = TRUE) +
  annotate("text", x=2, y=-0.05, label=paste("italic(P) == ", format.pval(lof_p, digits = 2)), parse = TRUE) +
  scale_x_discrete(labels = c('All variants', 'Variants in \n LoF intolerant genes'))

### Multi-tissue
tissue_all_variants <- fread(paste0(DATA, "MTSplice_all_tissues.csv"), drop=1, sep = ',')
tissue_all_variants <- tissue_all_variants[abs(delta_logit_psi) > 0.05]
tissue_all_variants <- melt(tissue_all_variants, measure.vars = tissues)
all_variants_mtsplice <- summarySE(tissue_all_variants, 
                                   measurevar = 'value', 
                                   groupvars = c('variable', 'Proband'))
lof_variants_mtsplice <- summarySE(tissue_all_variants[(lof_gene)], 
                                   measurevar = 'value', 
                                   groupvars = c('variable', 'Proband'))
variants_mtsplice <- rbindlist(list(all_variants_mtsplice, 
                                    lof_variants_mtsplice), 
                               idcol = 'group')
variants_mtsplice[, group := as.factor(group)]
variants_mtsplice[, Proband := ifelse(Proband, 'Proband', 'Sibling')]
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

variants_mtsplice$tissue <- tissues2short[variants_mtsplice$variable]
variants_mtsplice <- dcast(variants_mtsplice, group + tissue ~ Proband, value.var=c("value"))
variants_mtsplice[, diff := Proband-Sibling]
variants_mtsplice[, brain := grepl('Brain', tissue)]
variants_mtsplice[, brain := ifelse(brain, "Brain", "Other tissues")]
variants_mtsplice[, outlier := ifelse(is_outlier(diff), tissue, ""), by=c('group', 'brain')]
variants_mtsplice[, group := ifelse(group==1, 'All variants', 'Variants in LoF intolerant genes')]

brain_noneBrain <- ggplot(variants_mtsplice, aes(brain, diff)) + 
  geom_boxplot() +
  theme_classic() +
  geom_text_repel(aes(label=outlier),
                  nudge_x=-0.2, nudge_y = -0.002, size=3) +
  facet_wrap(~group) +
  labs(x=element_blank(), y='Proband-Sibling mean difference')  +
  theme_classic(base_size = 15)
library(ggpval)
brain_noneBrain <- add_pval(brain_noneBrain, fold_change = T,
                            pairs = list(c(1,2)),
                            test = 't.test')

large <- plot_grid(mmsplice, brain_noneBrain, labels = c("A", "B"), nrow = 1, rel_widths = c(1.3, 2))
ggsave("./data/Figure7.pdf", width = 11, height = 5.5)


### MMSplice < 0.05, to see whether MTSplice capture variants otherwise missed by MMSplice
tissue_all_variants <- fread(paste0(DATA, "MTSplice_all_tissues.csv"), drop=1, sep = ',')
tissue_all_variants <- tissue_all_variants[abs(delta_logit_psi) < 0.05]

pval_mmsplice <- wilcox.test(tissue_all_variants[Proband==T, delta_logit_psi],
                     tissue_all_variants[Proband==F, delta_logit_psi],
                     alternative='less')$p.value

pval_mtsplice <- wilcox.test(tissue_all_variants[Proband==T, "Frontal Cortex - Brain"][[1]],
                             tissue_all_variants[Proband==F, "Frontal Cortex - Brain"][[1]],
                             alternative='less')$p.value
mmsplice <- summarySE(tissue_all_variants, 
                      measurevar = 'delta_logit_psi', 
                      groupvars = c('Proband'))
setnames(mmsplice, 'delta_logit_psi', 'effect')
mtsplice <- summarySE(tissue_all_variants, 
                      measurevar = 'Frontal Cortex - Brain', 
                      groupvars = c('Proband'))
setnames(mtsplice, 'Frontal Cortex - Brain', 'effect')
df <- rbindlist(list(mmsplice, mtsplice), idcol = c('group'))

df[, group := as.factor(group)]
df[, Proband := ifelse(Proband, 'Proband', 'Sibling')]

# Compare MMSplice and MTSplice
comp <- ggplot(df, aes(group, effect, color=Proband)) +
  geom_pointrange(aes(ymin=effect-ci, ymax=effect+ci), 
                  position = position_dodge(width = 0.3)) +
  scale_colour_manual(values=cbPalette) +
  labs(x=element_blank(), y=expression('Predicted '*Delta*'logit'*Psi)) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(), legend.position = 'bottom') +
  annotate("text", x=1, y=0.003, label=paste("italic(P) == ", format.pval(pval_mmsplice, digits = 2)), parse = TRUE) +
  annotate("text", x=2, y=0.003, label=paste("italic(P) == ", format.pval(pval_mtsplice, digits = 2)), parse = TRUE) +
  scale_x_discrete(labels = c('MMSplice', 'MTSplice \n Frontal Cortex'))

# Compare Brain vs non-brain
tissue_all_variants <- melt(tissue_all_variants, measure.vars = tissues)
effects_mtsplice <- summarySE(tissue_all_variants, 
          measurevar = 'value', 
          groupvars = c('variable', 'Proband'))
effects_mtsplice <- as.data.table(effects_mtsplice)
effects_mtsplice[, Proband := ifelse(Proband, 'Proband', 'Sibling')]
effects_mtsplice$tissue <- tissues2short[effects_mtsplice$variable]
effects_mtsplice <- dcast(effects_mtsplice, tissue ~ Proband, value.var=c("value"))
effects_mtsplice[, diff := Proband-Sibling]
effects_mtsplice[, brain := grepl('Brain', tissue)]
effects_mtsplice[, brain := ifelse(brain, "Brain", "Other tissues")]
effects_mtsplice[, outlier := ifelse(is_outlier(diff), tissue, ""), by=c('brain')]

mtsplice_small <- ggplot(effects_mtsplice, aes(brain, diff)) + 
  geom_boxplot() +
  theme_classic() +
  geom_text_repel(aes(label=outlier),
                  nudge_x=-0.2, nudge_y = -0.001, size=3) +
  labs(x=element_blank(), y='Proband-Sibling mean difference')  +
  theme_classic(base_size = 15)

mtsplice_small <- add_pval(mtsplice_small, fold_change = F,
         pairs = list(c(1,2)),
         test = 't.test')

small <- plot_grid(comp, mtsplice_small, labels = c("C", "D"), nrow = 1, rel_widths = c(1.3, 2))

ggsave("./data/Sup_MMSplice_small.pdf", width = 11, height = 5.5)
