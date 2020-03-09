# (featureTable %in% mf_5sp_allSamples_forAPG.txt) # Make sure samples in feature table match mapping file?
# (species !=“not”)  # First, get rid of all non-amphibian samples
# (orig_contam==0)  # get rid of amphibians that were contaminated upon arrival
# (Bd_exposure==“Bd-exposed”) # Get ONLY exposed individuals (no controls)
# (time<=5) # Get pre-exposured time points (exposure occurred between time points 5 and 6)
# —> I would also recommend “combing” all time points by aggregating by indivID at this point, if you don’t want to violate independence assumptions.

# qiime feature-table filter-samples --i-table featureTable_count_wtaxa_filtered.qza --o-filtered-table non-rarefied.qza --m-metadata-file mf_5sp_allSamples_forAPG.txt --p-where "species != 'not' and orig_contam == '0' and Bd_exposure == 'Bd-exposed' and time < 5"
# qiime feature-table summarize --i-table non-rarefied.qza --o-visualization non-rarefied.qzv --m-sample-metadata-file mf_5sp_allSamples_forAPG.txt

# qiime feature-table filter-samples --i-table featureTable_r10k_wtaxa_filtered.qza --o-filtered-table rarefied.qza --m-metadata-file mf_5sp_allSamples_forAPG.txt --p-where "species != 'not' and orig_contam == '0' and Bd_exposure == 'Bd-exposed' and time < 5"
qiime feature-table summarize --i-table rarefied.qza --o-visualization rarefied.qzv --m-sample-metadata-file mf_5sp_allSamples_forAPG.txt

# Make stratified sampling version
# qiime feature-table filter-samples --i-table featureTable_count_wtaxa_filtered.qza --o-filtered-table stratified-nonrare.qza --m-metadata-file mf_stratsamp.txt 
# qiime feature-table summarize --i-table stratified-nonrare.qza --o-visualization stratified-nonrare.qzv --m-sample-metadata-file mf_stratsamp.txt


# Collapse samples by sum
# qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_summed.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table non-rarefied_summed.qza --o-visualization non-rarefied_summed.qzv --m-sample-metadata-file collapsed_metadata.txt

# qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_summed.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table rarefied_summed.qza --o-visualization rarefied_summed.qzv --m-sample-metadata-file collapsed_metadata.txt

# Collapse samples by median
# qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_median.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table non-rarefied_median.qza --o-visualization non-rarefied_median.qzv --m-sample-metadata-file collapsed_metadata.txt

# qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_median.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table rarefied_median.qza --o-visualization rarefied_median.qzv --m-sample-metadata-file collapsed_metadata.txt

# Collapse samples by mean
# qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_mean.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table non-rarefied_mean.qza --o-visualization non-rarefied_mean.qzv --m-sample-metadata-file collapsed_metadata.txt
# 
# qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_mean.qza \
# --p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
# --p-mode 'sum' 
# qiime feature-table summarize --i-table rarefied_mean.qza --o-visualization rarefied_mean.qzv --m-sample-metadata-file collapsed_metadata.txt


# Note to self: export predictions and feature_importance to plot in R manually.

# Non-rarefied, stratified
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.1/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.2/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.3/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.4/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.5/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.6/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir stratified-nonrare_0.7/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7
 
# qiime metadata tabulate \
#   --m-input-file stratified-nonrare/predictions.qza \
#   --o-visualization stratified-nonrare/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file stratified-nonrare/feature_importance.qza \
#   --o-visualization stratified-nonrare/feature_importance.qzv

# Non-rarefied, collapsed, sum
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.1/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.2/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.3/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.4/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.5/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.6/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir non-rarefied_summed_0.7/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

# qiime metadata tabulate \
#   --m-input-file non-rarefied_summed/predictions.qza \
#   --o-visualization non-rarefied_summed/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied_summed/feature_importance.qza \
#   --o-visualization non-rarefied_summed/feature_importance.qzv
  
  # Non-rarefied, collapsed, median
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.1/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.2/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.3/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.4/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.5/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.6/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir non-rarefied_median_0.7/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

# qiime metadata tabulate \
#   --m-input-file non-rarefied_median/predictions.qza \
#   --o-visualization non-rarefied_median/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied_median/feature_importance.qza \
#   --o-visualization non-rarefied_median/feature_importance.qzv
  
    
  # Non-rarefied, collapsed, mean
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.1/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.2/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.3/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.4/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.5/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.6/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir non-rarefied_mean_0.7/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

# qiime metadata tabulate \
#   --m-input-file non-rarefied_mean/predictions.qza \
#   --o-visualization non-rarefied_mean/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied_mean/feature_importance.qza \
#   --o-visualization non-rarefied_mean/feature_importance.qzv
  
# Note to self: export predictions and feature_importance to plot in R manually.
# Compare with non-collapsed prediction accuracy






  
#### Non-rarefied,binary, non-collapsed
# qiime sample-classifier classify-samples --output-dir non-rarefied-binary/ --i-table non-rarefied.qza \
#   --m-metadata-file mf_5sp_allSamples_forAPG_binary.txt --m-metadata-column future_max_bd_load_binary --p-n-estimators 1000 \
#   --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file non-rarefied-binary/predictions.qza \
#   --o-visualization non-rarefied-binary/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied-binary/feature_importance.qza \
#   --o-visualization non-rarefied-binary/feature_importance.qzv
# 

##### Non-rarefied,binary, collapsed
# qiime sample-classifier classify-samples --output-dir non-rarefied-binary_grouped/ --i-table non-rarefied_grouped.qza \
#   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load_binary --p-n-estimators 1000 \
#   --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file non-rarefied-binary_grouped/predictions.qza \
#   --o-visualization non-rarefied-binary_grouped/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied-binary_grouped/feature_importance.qza \
#   --o-visualization non-rarefied-binary_grouped/feature_importance.qzv



#### DOESn'T WORK; because feature table is not frequency or composition, it's a weird p-value thing.
# Community-level, collapsed
# biom convert -i community_level_traits.txt -o community_level_traits.biom --table-type="OTU table" --to-hdf5
# qiime tools import --type featureTable[Composition] --input-path community_level_traits.biom --output-path community_level_traits.qza
# 
# qiime sample-classifier regress-samples --output-dir community_level/ --i-table non-community_level_traits.qza \
#    --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
#    --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file community_level/predictions.qza \
#   --o-visualization community_level/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file community_level/feature_importance.qza \
#   --o-visualization community_level/feature_importance.qzv
#   
##### Note to self: export predictions and feature_importance to plot in R manually.
##### Compare with non-collapsed prediction accuracy




# Make a presence/absence table too-- 
# qiime feature-table presence-absence --i-table non-rarefied.qza \
# --o-presence-absence-table non-rarefied_PA.qza
# # Now, stratify
# qiime feature-table filter-samples --i-table non-rarefied_PA.qza --o-filtered-table stratified-PA.qza --m-metadata-file mf_stratsamp.txt 

# Manually made prevalence table
# biom convert -i OTUTable_prevalence.txt -o OTUTable_prevalence.biom --table-type "OTU table" --to-hdf5
# qiime tools import --i-table OTUTable_prevalence.biom




# Non-rarefied, non-collapsed
# qiime sample-classifier regress-samples --output-dir non-rarefied/ --i-table non-rarefied.qza \
#    --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
#    --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file non-rarefied/predictions.qza \
#   --o-visualization non-rarefied/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file non-rarefied/feature_importance.qza \
#   --o-visualization non-rarefied/feature_importance.qzv
#   


  
# Non-rarefied, PA, stratified--- to working because you need frequency
# qiime sample-classifier regress-samples --output-dir stratified-PA/ --i-table stratified-PA.qza \
#    --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
#    --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file stratified-PA/predictions.qza \
#   --o-visualization stratified-PA/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file stratified-PA/feature_importance.qza \
#   --o-visualization stratified-PA/feature_importance.qzv

  
  

##### Rarefied
# qiime sample-classifier regress-samples --output-dir rarefied/ --i-table rarefied.qza \
#    --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
#    --p-optimize-feature-selection --p-parameter-tuning
# qiime metadata tabulate \
#   --m-input-file rarefied/predictions.qza \
#   --o-visualization rarefied/predictions.qzv
# qiime metadata tabulate \
#   --m-input-file rarefied/feature_importance.qza \
#   --o-visualization rarefied/feature_importance.qzv
