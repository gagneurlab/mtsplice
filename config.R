library(data.table)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
theme_update(plot.title = element_text(hjust = 0.5))
library(magrittr)
library("rhdf5")

DATA <- './data/'
GTEX_DATA <- './data/processed/gtex/savedObjects'
GTEX_FRASER_PATH <- './data/processed/gtex/'
REPORT <- './reports/'

# function to get tissue name from snakemake name
get_tissue <- function(filen){
  splits <- strsplit(filen, '/')[[1]]
  tissue <- splits[length(splits) - 1]
  stopifnot(grepl('^gtex', tissue))
  tissue
}

# here I convert ascot tissue names to gtex tissue names, except the two eye tissues
tissues <- c('Retina - Eye', 'RPE/Choroid/Sclera - Eye', 'Subcutaneous - Adipose',
             'Visceral (Omentum) - Adipose', 'Adrenal Gland', 'Aorta - Artery',
             'Coronary - Artery', 'Tibial - Artery', 'Bladder', 'Amygdala - Brain',
             'Anterior cingulate - Brain', 'Caudate nucleus - Brain',
             'Cerebellar Hemisphere - Brain', 'Cerebellum - Brain', 'Cortex - Brain',
             'Frontal Cortex - Brain', 'Hippocampus - Brain', 'Hypothalamus - Brain',
             'Nucleus accumbens - Brain', 'Putamen - Brain',
             'Spinal cord (C1) - Brain', 'Substantia nigra - Brain',
             'Mammary Tissue - Breast', 'EBV-xform lymphocytes - Cells',
             'Leukemia (CML) - Cells', 'Xform. fibroblasts - Cells',
             'Ectocervix - Cervix', 'Endocervix - Cervix', 'Sigmoid - Colon',
             'Transverse - Colon', 'Gastroesoph. Junc. - Esophagus',
             'Mucosa - Esophagus', 'Muscularis - Esophagus', 'Fallopian Tube',
             'Atrial Appendage - Heart', 'Left Ventricle - Heart', 'Cortex - Kidney',
             'Liver', 'Lung', 'Minor Salivary Gland', 'Skeletal - Muscle',
             'Tibial - Nerve', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate',
             'Not Sun Exposed - Skin', 'Sun Exposed (Lower leg) - Skin',
             'Ileum - Small Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid',
             'Uterus', 'Vagina', 'Whole Blood')

# This is in GTEx tissue name format, some names are shorten
tissues_gtex <- c(
             'Adipose - Subcutaneous',
             'Adipose - Visceral (Omentum)',
             'Adrenal Gland',
             'Artery - Aorta',
             'Artery - Coronary',
             'Artery - Tibial',
             'Bladder',
             'Brain - Amygdala',
             'Brain - Anterior cingulate cortex',
             'Brain - Caudate',
             'Brain - Cerebellar Hemisphere',
             'Brain - Cerebellum',
             'Brain - Cortex',
             'Brain - Frontal Cortex',
             'Brain - Hippocampus',
             'Brain - Hypothalamus',
             'Brain - Nucleus accumbens',
             'Brain - Putamen ',
             'Brain - Spinal cord',
             'Brain - Substantia nigra',
             'Breast - Mammary Tissue',
             'EBV-lymphocytes',
             'CML',
             'fibroblasts',
             'Endocervix',
             'Colon - Sigmoid',
             'Colon - Transverse',
             'Esophagus - G. Junction',
             'Esophagus - Mucosa',
             'Esophagus - Muscularis',
             'Fallopian Tube',
             'Heart - Atrial Appendage',
             'Heart - Left Ventricle',
             'Kidney - Cortex',
             'Liver',
             'Lung',
             'Minor Salivary Gland',
             'Muscle Skeletal',
             'Nerve - Tibial',
             'Ovary',
             'Pancreas',
             'Pituitary',
             'Prostate',
             'Skin - Not Sun Exposed',
             'Skin - Sun',
             'Small Intestine - Terminal Ileum',
             'Spleen',
             'Stomach',
             'Testis',
             'Thyroid',
             'Uterus',
             'Vagina',
             'Whole Blood')

# This is in GTEx tissue name format, some names are shorten
tissues_short <- c(
  'Retina - Eye', 
  'RPE/Choroid/Sclera - Eye',
  'Adipose - Subcutaneous',
  'Adipose - Visceral (Omentum)',
  'Adrenal Gland',
  'Artery - Aorta',
  'Artery - Coronary',
  'Artery - Tibial',
  'Bladder',
  'Brain - Amygdala',
  'Brain - Anterior cingulate cortex',
  'Brain - Caudate',
  'Brain - Cerebellar Hemisphere',
  'Brain - Cerebellum',
  'Brain - Cortex',
  'Brain - Frontal Cortex',
  'Brain - Hippocampus',
  'Brain - Hypothalamus',
  'Brain - Nucleus accumbens',
  'Brain - Putamen ',
  'Brain - Spinal cord',
  'Brain - Substantia nigra',
  'Breast - Mammary Tissue',
  'EBV-lymphocytes',
  'CML',
  'fibroblasts',
  'Ectocervix',
  'Endocervix',
  'Colon - Sigmoid',
  'Colon - Transverse',
  'Esophagus - G. Junction',
  'Esophagus - Mucosa',
  'Esophagus - Muscularis',
  'Fallopian Tube',
  'Heart - Atrial Appendage',
  'Heart - Left Ventricle',
  'Kidney - Cortex',
  'Liver',
  'Lung',
  'Minor Salivary Gland',
  'Muscle Skeletal',
  'Nerve - Tibial',
  'Ovary',
  'Pancreas',
  'Pituitary',
  'Prostate',
  'Skin - Not Sun Exposed',
  'Skin - Sun',
  'Small Intestine - Terminal Ileum',
  'Spleen',
  'Stomach',
  'Testis',
  'Thyroid',
  'Uterus',
  'Vagina',
  'Whole Blood')

tissues2short <- tissues_short
names(tissues2short) <- tissues

#
mse <- function(x, y){
  mean((x - y) ^ 2, na.rm=T)
}

rmse <- function(x, y){
  round(sqrt(mse(x,y)), 4)
}

nacor <- function(x, y, ...){
  cor(x, y, use="pairwise.complete.obs", ...)
}


#------------
# Format
format_R <- function(clt){
  clt <- round(clt, 2)
  bquote(italic(R) == .(format.pval(clt)))
}

format_mse <- function(clt){
  clt <- round(clt, 3)
  paste0("RMSE", "=", format.pval(clt))
}



summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

clip <- function(x, lower, upper){
  x[which(x < lower)] <- lower
  x[which(x > upper)] <- upper
  x
}

logit <- function(x){
  x <- clip(x, 1e-5, 1-1e-5)
  x / (1-x)
}

sigmoid <- function(x){
  1 / (1 + exp(-x))
}


format_pval <- function(pval){
  if (is.character(pval)){
    # pval contains fold change
    pval <- strsplit(pval, ' ')[[1]]
    fc <- pval[2]
    pval <- as.numeric(pval[1])
  }else{
    fc <- ""
  }
  pval <- format.pval(pval, digits = 2)
  if (grepl("<", pval)){
    pval <- gsub("< ?", "", pval)
    pval <- bquote(italic(P) < .(paste(pval, fc)))
  }else{
    pval <- bquote(italic(P) == .(paste(pval, fc)))
  }
  pval
}

