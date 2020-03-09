# (featureTable %in% mf_5sp_allSamples_forAPG.txt) # Make sure samples in feature table match mapping file?
# (species !=“not”)  # First, get rid of all non-amphibian samples
# (orig_contam==0)  # get rid of amphibians that were contaminated upon arrival
# (Bd_exposure==“Bd-exposed”) # Get ONLY exposed individuals (no controls)
# (time<=5) # Get pre-exposured time points (exposure occurred between time points 5 and 6)
# —> I would also recommend “combing” all time points by aggregating by indivID at this point, if you don’t want to violate independence assumptions.

qiime feature-table filter-samples --i-table featureTable_count_wtaxa_filtered.qza --o-filtered-table non-rarefied.qza --m-metadata-file mf_5sp_allSamples_forAPG.txt --p-where "species != 'not' and orig_contam == '0' and Bd_exposure == 'Bd-exposed' and time < 5"
qiime feature-table summarize --i-table non-rarefied.qza --o-visualization non-rarefied.qzv --m-sample-metadata-file mf_5sp_allSamples_forAPG.txt

qiime feature-table filter-samples --i-table featureTable_r10k_wtaxa_filtered.qza --o-filtered-table rarefied.qza --m-metadata-file mf_5sp_allSamples_forAPG.txt --p-where "species != 'not' and orig_contam == '0' and Bd_exposure == 'Bd-exposed' and time < 5"
qiime feature-table summarize --i-table rarefied.qza --o-visualization rarefied.qzv --m-sample-metadata-file mf_5sp_allSamples_forAPG.txt

####### Collapsing and Stratifying #########

# Make stratified sampling version-nonrare
qiime feature-table filter-samples --i-table featureTable_count_wtaxa_filtered.qza --o-filtered-table stratified-nonrare.qza --m-metadata-file mf_stratsamp.txt 
qiime feature-table summarize --i-table stratified-nonrare.qza --o-visualization stratified-nonrare.qzv --m-sample-metadata-file mf_stratsamp.txt

# Make stratified sampling version- rare
qiime feature-table filter-samples --i-table featureTable_r10k_wtaxa_filtered.qza --o-filtered-table stratified-rare.qza --m-metadata-file mf_stratsamp.txt 
qiime feature-table summarize --i-table stratified-rare.qza --o-visualization stratified-rare.qzv --m-sample-metadata-file mf_stratsamp.txt


# Collapse samples by sum
qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_summed.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table non-rarefied_summed.qza --o-visualization non-rarefied_summed.qzv --m-sample-metadata-file collapsed_metadata.txt

qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_summed.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table rarefied_summed.qza --o-visualization rarefied_summed.qzv --m-sample-metadata-file collapsed_metadata.txt

# Collapse samples by median
qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_median.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table non-rarefied_median.qza --o-visualization non-rarefied_median.qzv --m-sample-metadata-file collapsed_metadata.txt

qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_median.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table rarefied_median.qza --o-visualization rarefied_median.qzv --m-sample-metadata-file collapsed_metadata.txt

# Collapse samples by mean
qiime feature-table group --i-table non-rarefied.qza --o-grouped-table non-rarefied_mean.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table non-rarefied_mean.qza --o-visualization non-rarefied_mean.qzv --m-sample-metadata-file collapsed_metadata.txt

qiime feature-table group --i-table rarefied.qza --o-grouped-table rarefied_mean.qza \
--p-axis 'sample' --m-metadata-file mf_5sp_allSamples_forAPG.txt --m-metadata-column indivID \
--p-mode 'sum' 
qiime feature-table summarize --i-table rarefied_mean.qza --o-visualization rarefied_mean.qzv --m-sample-metadata-file collapsed_metadata.txt

####### Non-rarefied #########

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

####### Rarefied #########

# Rarefied, stratified
qiime sample-classifier regress-samples --output-dir stratified-rare_0.1/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir stratified-rare_0.2/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir stratified-rare_0.3/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir stratified-rare_0.4/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir stratified-rare_0.5/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir stratified-rare_0.6/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir stratified-rare_0.7/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

# Rarefied, collapsed, sum
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.1/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.2/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.3/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.4/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.5/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.6/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir rarefied_summed_0.7/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7


  # Rarefied, collapsed, median
qiime sample-classifier regress-samples --output-dir rarefied_median_0.1/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir rarefied_median_0.2/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir rarefied_median_0.3/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir rarefied_median_0.4/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir rarefied_median_0.5/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir rarefied_median_0.6/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir rarefied_median_0.7/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

  
    
  # Rarefied, collapsed, mean
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.1/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.2/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.3/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.4/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.5/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.6/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier regress-samples --output-dir rarefied_mean_0.7/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column future_max_bd_load --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7


######## Exporting ##########

# Non-rarefied, stratified
mkdir ALL_stratified-nonrare
qiime tools export --input-path stratified-nonrare_0.1/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.1
qiime tools export --input-path stratified-nonrare_0.1/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.1

qiime tools export --input-path stratified-nonrare_0.2/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.2
qiime tools export --input-path stratified-nonrare_0.2/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.2

qiime tools export --input-path stratified-nonrare_0.3/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.3
qiime tools export --input-path stratified-nonrare_0.3/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.3
   
qiime tools export --input-path stratified-nonrare_0.4/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.4
qiime tools export --input-path stratified-nonrare_0.4/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.4

qiime tools export --input-path stratified-nonrare_0.5/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.5
qiime tools export --input-path stratified-nonrare_0.5/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.5

qiime tools export --input-path stratified-nonrare_0.6/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.6
qiime tools export --input-path stratified-nonrare_0.6/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.6

qiime tools export --input-path stratified-nonrare_0.7/feature_importance.qza --output-path ALL_stratified-nonrare/features_0.7
qiime tools export --input-path stratified-nonrare_0.7/predictions.qza --output-path ALL_stratified-nonrare/predictions_0.7


mkdir ALL_stratified-rare
qiime tools export --input-path stratified-rare_0.1/feature_importance.qza --output-path ALL_stratified-rare/features_0.1
qiime tools export --input-path stratified-rare_0.1/predictions.qza --output-path ALL_stratified-rare/predictions_0.1

qiime tools export --input-path stratified-rare_0.2/feature_importance.qza --output-path ALL_stratified-rare/features_0.2
qiime tools export --input-path stratified-rare_0.2/predictions.qza --output-path ALL_stratified-rare/predictions_0.2

qiime tools export --input-path stratified-rare_0.3/feature_importance.qza --output-path ALL_stratified-rare/features_0.3
qiime tools export --input-path stratified-rare_0.3/predictions.qza --output-path ALL_stratified-rare/predictions_0.3
   
qiime tools export --input-path stratified-rare_0.4/feature_importance.qza --output-path ALL_stratified-rare/features_0.4
qiime tools export --input-path stratified-rare_0.4/predictions.qza --output-path ALL_stratified-rare/predictions_0.4

qiime tools export --input-path stratified-rare_0.5/feature_importance.qza --output-path ALL_stratified-rare/features_0.5
qiime tools export --input-path stratified-rare_0.5/predictions.qza --output-path ALL_stratified-rare/predictions_0.5

qiime tools export --input-path stratified-rare_0.6/feature_importance.qza --output-path ALL_stratified-rare/features_0.6
qiime tools export --input-path stratified-rare_0.6/predictions.qza --output-path ALL_stratified-rare/predictions_0.6

qiime tools export --input-path stratified-rare_0.7/feature_importance.qza --output-path ALL_stratified-rare/features_0.7
qiime tools export --input-path stratified-rare_0.7/predictions.qza --output-path ALL_stratified-rare/predictions_0.7


# Non-rarefied, collapsed, sum
mkdir ALL_non-rarefied_sum

qiime tools export --input-path non-rarefied_summed_0.1/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.1
qiime tools export --input-path non-rarefied_summed_0.1/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.1

qiime tools export --input-path non-rarefied_summed_0.2/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.2
qiime tools export --input-path non-rarefied_summed_0.2/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.2

qiime tools export --input-path non-rarefied_summed_0.3/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.3
qiime tools export --input-path non-rarefied_summed_0.3/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.3
   
qiime tools export --input-path non-rarefied_summed_0.4/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.4
qiime tools export --input-path non-rarefied_summed_0.4/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.4

qiime tools export --input-path non-rarefied_summed_0.5/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.5
qiime tools export --input-path non-rarefied_summed_0.5/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.5

qiime tools export --input-path non-rarefied_summed_0.6/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.6
qiime tools export --input-path non-rarefied_summed_0.6/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.6

qiime tools export --input-path non-rarefied_summed_0.7/feature_importance.qza --output-path ALL_non-rarefied_sum/features_0.7
qiime tools export --input-path non-rarefied_summed_0.7/predictions.qza --output-path ALL_non-rarefied_sum/predictions_0.7

# Rarefied, collapsed, sum
mkdir ALL_rarefied_sum

qiime tools export --input-path rarefied_summed_0.1/feature_importance.qza --output-path ALL_rarefied_sum/features_0.1
qiime tools export --input-path rarefied_summed_0.1/predictions.qza --output-path ALL_rarefied_sum/predictions_0.1

qiime tools export --input-path rarefied_summed_0.2/feature_importance.qza --output-path ALL_rarefied_sum/features_0.2
qiime tools export --input-path rarefied_summed_0.2/predictions.qza --output-path ALL_rarefied_sum/predictions_0.2

qiime tools export --input-path rarefied_summed_0.3/feature_importance.qza --output-path ALL_rarefied_sum/features_0.3
qiime tools export --input-path rarefied_summed_0.3/predictions.qza --output-path ALL_rarefied_sum/predictions_0.3
   
qiime tools export --input-path rarefied_summed_0.4/feature_importance.qza --output-path ALL_rarefied_sum/features_0.4
qiime tools export --input-path rarefied_summed_0.4/predictions.qza --output-path ALL_rarefied_sum/predictions_0.4

qiime tools export --input-path rarefied_summed_0.5/feature_importance.qza --output-path ALL_rarefied_sum/features_0.5
qiime tools export --input-path rarefied_summed_0.5/predictions.qza --output-path ALL_rarefied_sum/predictions_0.5

qiime tools export --input-path rarefied_summed_0.6/feature_importance.qza --output-path ALL_rarefied_sum/features_0.6
qiime tools export --input-path rarefied_summed_0.6/predictions.qza --output-path ALL_rarefied_sum/predictions_0.6

qiime tools export --input-path rarefied_summed_0.7/feature_importance.qza --output-path ALL_rarefied_sum/features_0.7
qiime tools export --input-path rarefied_summed_0.7/predictions.qza --output-path ALL_rarefied_sum/predictions_0.7


  
# Non-rarefied, collapsed, median
mkdir ALL_non-rarefied_median

qiime tools export --input-path non-rarefied_median_0.1/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.1
qiime tools export --input-path non-rarefied_median_0.1/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.1

qiime tools export --input-path non-rarefied_median_0.2/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.2
qiime tools export --input-path non-rarefied_median_0.2/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.2

qiime tools export --input-path non-rarefied_median_0.3/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.3
qiime tools export --input-path non-rarefied_median_0.3/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.3
   
qiime tools export --input-path non-rarefied_median_0.4/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.4
qiime tools export --input-path non-rarefied_median_0.4/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.4

qiime tools export --input-path non-rarefied_median_0.5/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.5
qiime tools export --input-path non-rarefied_median_0.5/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.5

qiime tools export --input-path non-rarefied_median_0.6/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.6
qiime tools export --input-path non-rarefied_median_0.6/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.6

qiime tools export --input-path non-rarefied_median_0.7/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.7
qiime tools export --input-path non-rarefied_median_0.7/predictions.qza --output-path ALL_non-rarefied_median/predictions_0.7

# Rarefied, collapsed, median
mkdir ALL_rarefied_median

qiime tools export --input-path rarefied_median_0.1/feature_importance.qza --output-path ALL_rarefied_median/features_0.1
qiime tools export --input-path rarefied_median_0.1/predictions.qza --output-path ALL_rarefied_median/predictions_0.1

qiime tools export --input-path rarefied_median_0.2/feature_importance.qza --output-path ALL_rarefied_median/features_0.2
qiime tools export --input-path rarefied_median_0.2/predictions.qza --output-path ALL_rarefied_median/predictions_0.2

qiime tools export --input-path rarefied_median_0.3/feature_importance.qza --output-path ALL_rarefied_median/features_0.3
qiime tools export --input-path rarefied_median_0.3/predictions.qza --output-path ALL_rarefied_median/predictions_0.3
   
qiime tools export --input-path rarefied_median_0.4/feature_importance.qza --output-path ALL_rarefied_median/features_0.4
qiime tools export --input-path rarefied_median_0.4/predictions.qza --output-path ALL_rarefied_median/predictions_0.4

qiime tools export --input-path rarefied_median_0.5/feature_importance.qza --output-path ALL_rarefied_median/features_0.5
qiime tools export --input-path rarefied_median_0.5/predictions.qza --output-path ALL_rarefied_median/predictions_0.5

qiime tools export --input-path rarefied_median_0.6/feature_importance.qza --output-path ALL_rarefied_median/features_0.6
qiime tools export --input-path rarefied_median_0.6/predictions.qza --output-path ALL_rarefied_median/predictions_0.6

qiime tools export --input-path rarefied_median_0.7/feature_importance.qza --output-path ALL_rarefied_median/features_0.7
qiime tools export --input-path rarefied_median_0.7/predictions.qza --output-path ALL_rarefied_median/predictions_0.7

# Non-rarefied, collapsed, mean
mkdir ALL_non-rarefied_mean
qiime tools export --input-path non-rarefied_mean_0.1/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.1
qiime tools export --input-path non-rarefied_mean_0.1/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.1

qiime tools export --input-path non-rarefied_mean_0.2/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.2
qiime tools export --input-path non-rarefied_mean_0.2/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.2

qiime tools export --input-path non-rarefied_mean_0.3/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.3
qiime tools export --input-path non-rarefied_mean_0.3/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.3
   
qiime tools export --input-path non-rarefied_mean_0.4/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.4
qiime tools export --input-path non-rarefied_mean_0.4/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.4

qiime tools export --input-path non-rarefied_mean_0.5/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.5
qiime tools export --input-path non-rarefied_mean_0.5/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.5

qiime tools export --input-path non-rarefied_mean_0.6/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.6
qiime tools export --input-path non-rarefied_mean_0.6/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.6

qiime tools export --input-path non-rarefied_mean_0.7/feature_importance.qza --output-path ALL_non-rarefied_mean/features_0.7
qiime tools export --input-path non-rarefied_mean_0.7/predictions.qza --output-path ALL_non-rarefied_mean/predictions_0.7


# Rarefied, collapsed, mean
mkdir ALL_rarefied_mean
qiime tools export --input-path rarefied_mean_0.1/feature_importance.qza --output-path ALL_rarefied_mean/features_0.1
qiime tools export --input-path rarefied_mean_0.1/predictions.qza --output-path ALL_rarefied_mean/predictions_0.1

qiime tools export --input-path rarefied_mean_0.2/feature_importance.qza --output-path ALL_rarefied_mean/features_0.2
qiime tools export --input-path non-rarefied_mean_0.2/predictions.qza --output-path ALL_rarefied_mean/predictions_0.2

qiime tools export --input-path rarefied_mean_0.3/feature_importance.qza --output-path ALL_rarefied_mean/features_0.3
qiime tools export --input-path rarefied_mean_0.3/predictions.qza --output-path ALL_rarefied_mean/predictions_0.3
   
qiime tools export --input-path rarefied_mean_0.4/feature_importance.qza --output-path ALL_rarefied_mean/features_0.4
qiime tools export --input-path rarefied_mean_0.4/predictions.qza --output-path ALL_rarefied_mean/predictions_0.4

qiime tools export --input-path rarefied_mean_0.5/feature_importance.qza --output-path ALL_rarefied_mean/features_0.5
qiime tools export --input-path rarefied_mean_0.5/predictions.qza --output-path ALL_rarefied_mean/predictions_0.5

qiime tools export --input-path rarefied_mean_0.6/feature_importance.qza --output-path ALL_rarefied_mean/features_0.6
qiime tools export --input-path rarefied_mean_0.6/predictions.qza --output-path ALL_rarefied_mean/predictions_0.6

qiime tools export --input-path rarefied_mean_0.7/feature_importance.qza --output-path ALL_rarefied_mean/features_0.7
qiime tools export --input-path rarefied_mean_0.7/predictions.qza --output-path ALL_rarefied_mean/predictions_0.7

