#!bin/bash RScript

library("tidyverse")
#### Running decontam on dataset ###
# Load files
otu <- read.delim("./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt"
                  , skip=1
                  , header=TRUE
                  , row.names = 1)
mf <- read.delim("./0_QIITA_files/intermediate_files/mapping_5Sp_fixed.txt"
                 , header=TRUE)

# Transpost otu table, and get rid of X. at beignning
sampleNames <- gsub("^X", "", colnames(otu))
colnames(otu) <- sampleNames

# Get taxa
taxa <- otu[,'taxonomy']
names(taxa) <- rownames(otu)
otu.notaxa <- otu[,-ncol(otu)]

# Calculate otu richness in each sample
otu.PA <- otu.notaxa
otu.PA[otu.PA>0] <- 1
sampleRich <- colSums(otu.PA)
# Correlate OTU richness with abundance of each OTU
corrOTU <- data.frame(otu = names(taxa),taxa=as.vector(taxa), corPvalue=rep(NA,length=length(taxa)), t=rep(NA,length=length(taxa)))
for ( o in 1:length(taxa) ) {
  corrOTU[o,'corPvalue'] <- cor.test(as.numeric(otu.notaxa[o,]), sampleRich)$p.value
  corrOTU[o,'t'] <- cor.test(as.numeric(otu.notaxa[o,]), sampleRich)$statistic
}

potential_contam <- corrOTU %>%
  mutate(sig = corPvalue < 0.05, dir = t < 0) %>%
  filter(sig & dir)

# Now, filter out OTUs that are found in "controls"- aka, water
controls <- mf %>%
  filter(sample_type %in% c("water","media"), X.SampleID %in% sampleNames) %>%
  select(X.SampleID) %>%
  pull()

# Get otus that are found in controls
control_otus <- otu.PA %>%
  select(as.character(controls)) %>%
  rowSums()
contam_otus_in_controls <- names(control_otus)[(which(control_otus>0))]

examine_contam <- potential_contam %>%
  filter(otu %in% contam_otus_in_controls) %>%
  select(otu) %>%
  pull()

dir.create("Possible_contaminants")
for ( o in as.character(examine_contam) ) {
  count_otu <- as.vector(otu[o,-ncol(otu)])
  tax <- as.character(otu[o,'taxonomy'])
  ggsave(file=paste0("Possible_contaminants/",tax,".pdf"), 
  data.frame(Richness=as.numeric(sampleRich), count=as.numeric(count_otu)) %>%
    ggplot(aes(x=Richness, y=count)) +
    geom_point() +
    ggtitle(paste0(tax), subtitle = o)+
    theme(title =element_text(size=5))
  )
}


