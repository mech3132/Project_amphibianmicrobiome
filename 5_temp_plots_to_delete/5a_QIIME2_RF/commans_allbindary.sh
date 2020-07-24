# qiime sample-classifier classify-samples --output-dir non-rarefied-binary/ --i-table non-rarefied.qza \
#   --m-metadata-file mf_5sp_allSamples_forAPG_binary.txt --m-metadata-column future_max_bd_load_binary --p-n-estimators 1000 \
#   --p-optimize-feature-selection --p-parameter-tuning

    ####### Non-rarefied #########

# Non-rarefied, stratified
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.1/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.2/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.3/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.4/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.5/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.6/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir stratified-nonrare_bin_0.7/ --i-table stratified-nonrare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7
 
# Non-rarefied, collapsed, sum
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.1/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.2/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.3/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.4/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.5/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.6/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir non-rarefied_summed_bin_0.7/ --i-table non-rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

  # Non-rarefied, collapsed, median
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.1/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.2/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.3/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.4/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.5/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.6/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir non-rarefied_median_bin_0.7/ --i-table non-rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7 
    
  # Non-rarefied, collapsed, mean
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.1/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.2/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.3/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.4/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.5/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.6/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir non-rarefied_mean_bin_0.7/ --i-table non-rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

####### Rarefied #########

# Rarefied, stratified
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.1/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.2/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.3/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.4/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.5/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.6/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir stratified-rare_bin_0.7/ --i-table stratified-rare.qza \
   --m-metadata-file mf_stratsamp.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

# Rarefied, collapsed, sum
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.1/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.2/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.3/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.4/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.5/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.6/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir rarefied_summed_bin_0.7/ --i-table rarefied_summed.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7


  # Rarefied, collapsed, median
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.1/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.2/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.3/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.4/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.5/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.6/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir rarefied_median_bin_0.7/ --i-table rarefied_median.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7

  
    
  # Rarefied, collapsed, mean
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.1/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.1
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.2/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.2
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.3/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.3
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.4/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.4
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.5/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.5
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.6/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.6
qiime sample-classifier classify-samples --output-dir rarefied_mean_bin_0.7/ --i-table rarefied_mean.qza \
   --m-metadata-file collapsed_metadata.txt --m-metadata-column presence_absence --p-n-estimators 1000 \
   --p-optimize-feature-selection --p-parameter-tuning --p-test-size 0.7



######## Exporting ##########

# Non-rarefied, stratified
mkdir ALL_stratified-nonrare_bin
qiime tools export --input-path stratified-nonrare_bin_0.1/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.1
qiime tools export --input-path stratified-nonrare_bin_0.1/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.1

qiime tools export --input-path stratified-nonrare_bin_0.2/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.2
qiime tools export --input-path stratified-nonrare_bin_0.2/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.2

qiime tools export --input-path stratified-nonrare_bin_0.3/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.3
qiime tools export --input-path stratified-nonrare_bin_0.3/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.3
   
qiime tools export --input-path stratified-nonrare_bin_0.4/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.4
qiime tools export --input-path stratified-nonrare_bin_0.4/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.4

qiime tools export --input-path stratified-nonrare_bin_0.5/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.5
qiime tools export --input-path stratified-nonrare_bin_0.5/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.5

qiime tools export --input-path stratified-nonrare_bin_0.6/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.6
qiime tools export --input-path stratified-nonrare_bin_0.6/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.6

qiime tools export --input-path stratified-nonrare_bin_0.7/feature_importance.qza --output-path ALL_stratified-nonrare_bin/features_0.7
qiime tools export --input-path stratified-nonrare_bin_0.7/predictions.qza --output-path ALL_stratified-nonrare_bin/predictions_0.7


mkdir ALL_stratified-rare_bin
qiime tools export --input-path stratified-rare_bin_0.1/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.1
qiime tools export --input-path stratified-rare_bin_0.1/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.1

qiime tools export --input-path stratified-rare_bin_0.2/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.2
qiime tools export --input-path stratified-rare_bin_0.2/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.2

qiime tools export --input-path stratified-rare_bin_0.3/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.3
qiime tools export --input-path stratified-rare_bin_0.3/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.3
   
qiime tools export --input-path stratified-rare_bin_0.4/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.4
qiime tools export --input-path stratified-rare_bin_0.4/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.4

qiime tools export --input-path stratified-rare_bin_0.5/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.5
qiime tools export --input-path stratified-rare_bin_0.5/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.5

qiime tools export --input-path stratified-rare_bin_0.6/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.6
qiime tools export --input-path stratified-rare_bin_0.6/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.6

qiime tools export --input-path stratified-rare_bin_0.7/feature_importance.qza --output-path ALL_stratified-rare_bin/features_0.7
qiime tools export --input-path stratified-rare_bin_0.7/predictions.qza --output-path ALL_stratified-rare_bin/predictions_0.7


# Non-rarefied, collapsed, sum
mkdir ALL_non-rarefied_sum_bin

qiime tools export --input-path non-rarefied_summed_bin_0.1/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.1
qiime tools export --input-path non-rarefied_summed_bin_0.1/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.1

qiime tools export --input-path non-rarefied_summed_bin_0.2/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.2
qiime tools export --input-path non-rarefied_summed_bin_0.2/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.2

qiime tools export --input-path non-rarefied_summed_bin_0.3/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.3
qiime tools export --input-path non-rarefied_summed_bin_0.3/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.3
   
qiime tools export --input-path non-rarefied_summed_bin_0.4/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.4
qiime tools export --input-path non-rarefied_summed_bin_0.4/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.4

qiime tools export --input-path non-rarefied_summed_bin_0.5/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.5
qiime tools export --input-path non-rarefied_summed_bin_0.5/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.5

qiime tools export --input-path non-rarefied_summed_bin_0.6/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.6
qiime tools export --input-path non-rarefied_summed_bin_0.6/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.6

qiime tools export --input-path non-rarefied_summed_bin_0.7/feature_importance.qza --output-path ALL_non-rarefied_sum_bin/features_0.7
qiime tools export --input-path non-rarefied_summed_bin_0.7/predictions.qza --output-path ALL_non-rarefied_sum_bin/predictions_0.7

# Rarefied, collapsed, sum
mkdir ALL_rarefied_sum_bin

qiime tools export --input-path rarefied_summed_bin_0.1/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.1
qiime tools export --input-path rarefied_summed_bin_0.1/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.1

qiime tools export --input-path rarefied_summed_bin_0.2/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.2
qiime tools export --input-path rarefied_summed_bin_0.2/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.2

qiime tools export --input-path rarefied_summed_bin_0.3/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.3
qiime tools export --input-path rarefied_summed_bin_0.3/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.3
   
qiime tools export --input-path rarefied_summed_bin_0.4/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.4
qiime tools export --input-path rarefied_summed_bin_0.4/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.4

qiime tools export --input-path rarefied_summed_bin_0.5/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.5
qiime tools export --input-path rarefied_summed_bin_0.5/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.5

qiime tools export --input-path rarefied_summed_bin_0.6/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.6
qiime tools export --input-path rarefied_summed_bin_0.6/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.6

qiime tools export --input-path rarefied_summed_bin_0.7/feature_importance.qza --output-path ALL_rarefied_sum_bin/features_0.7
qiime tools export --input-path rarefied_summed_bin_0.7/predictions.qza --output-path ALL_rarefied_sum_bin/predictions_0.7


  
# Non-rarefied, collapsed, median
mkdir ALL_non-rarefied_median_bin

qiime tools export --input-path non-rarefied_median_bin_0.1/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.1
qiime tools export --input-path non-rarefied_median_bin_0.1/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.1

qiime tools export --input-path non-rarefied_median_bin_0.2/feature_importance.qza --output-path ALL_non-rarefied_median/features_0.2
qiime tools export --input-path non-rarefied_median_bin_0.2/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.2

qiime tools export --input-path non-rarefied_median_bin_0.3/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.3
qiime tools export --input-path non-rarefied_median_bin_0.3/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.3
   
qiime tools export --input-path non-rarefied_median_bin_0.4/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.4
qiime tools export --input-path non-rarefied_median_bin_0.4/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.4

qiime tools export --input-path non-rarefied_median_bin_0.5/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.5
qiime tools export --input-path non-rarefied_median_bin_0.5/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.5

qiime tools export --input-path non-rarefied_median_bin_0.6/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.6
qiime tools export --input-path non-rarefied_median_bin_0.6/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.6

qiime tools export --input-path non-rarefied_median_bin_0.7/feature_importance.qza --output-path ALL_non-rarefied_median_bin/features_0.7
qiime tools export --input-path non-rarefied_median_bin_0.7/predictions.qza --output-path ALL_non-rarefied_median_bin/predictions_0.7

# Rarefied, collapsed, median
mkdir ALL_rarefied_median_bin

qiime tools export --input-path rarefied_median_bin_0.1/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.1
qiime tools export --input-path rarefied_median_bin_0.1/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.1

qiime tools export --input-path rarefied_median_bin_0.2/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.2
qiime tools export --input-path rarefied_median_bin_0.2/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.2

qiime tools export --input-path rarefied_median_bin_0.3/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.3
qiime tools export --input-path rarefied_median_bin_0.3/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.3
   
qiime tools export --input-path rarefied_median_bin_0.4/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.4
qiime tools export --input-path rarefied_median_bin_0.4/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.4

qiime tools export --input-path rarefied_median_bin_0.5/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.5
qiime tools export --input-path rarefied_median_bin_0.5/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.5

qiime tools export --input-path rarefied_median_bin_0.6/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.6
qiime tools export --input-path rarefied_median_bin_0.6/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.6

qiime tools export --input-path rarefied_median_bin_0.7/feature_importance.qza --output-path ALL_rarefied_median_bin/features_0.7
qiime tools export --input-path rarefied_median_bin_0.7/predictions.qza --output-path ALL_rarefied_median_bin/predictions_0.7

# Non-rarefied, collapsed, mean
mkdir ALL_non-rarefied_mean_bin
qiime tools export --input-path non-rarefied_mean_bin_0.1/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.1
qiime tools export --input-path non-rarefied_mean_bin_0.1/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.1

qiime tools export --input-path non-rarefied_mean_bin_0.2/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.2
qiime tools export --input-path non-rarefied_mean_bin_0.2/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.2

qiime tools export --input-path non-rarefied_mean_bin_0.3/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.3
qiime tools export --input-path non-rarefied_mean_bin_0.3/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.3
   
qiime tools export --input-path non-rarefied_mean_bin_0.4/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.4
qiime tools export --input-path non-rarefied_mean_bin_0.4/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.4

qiime tools export --input-path non-rarefied_mean_bin_0.5/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.5
qiime tools export --input-path non-rarefied_mean_bin_0.5/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.5

qiime tools export --input-path non-rarefied_mean_bin_0.6/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.6
qiime tools export --input-path non-rarefied_mean_bin_0.6/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.6

qiime tools export --input-path non-rarefied_mean_bin_0.7/feature_importance.qza --output-path ALL_non-rarefied_mean_bin/features_0.7
qiime tools export --input-path non-rarefied_mean_bin_0.7/predictions.qza --output-path ALL_non-rarefied_mean_bin/predictions_0.7


# Rarefied, collapsed, mean
mkdir ALL_rarefied_mean_bin
qiime tools export --input-path rarefied_mean_bin_0.1/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.1
qiime tools export --input-path rarefied_mean_bin_0.1/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.1

qiime tools export --input-path rarefied_mean_bin_0.2/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.2
qiime tools export --input-path non-rarefied_mean_bin_0.2/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.2

qiime tools export --input-path rarefied_mean_bin_0.3/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.3
qiime tools export --input-path rarefied_mean_bin_0.3/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.3
   
qiime tools export --input-path rarefied_mean_bin_0.4/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.4
qiime tools export --input-path rarefied_mean_bin_0.4/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.4

qiime tools export --input-path rarefied_mean_bin_0.5/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.5
qiime tools export --input-path rarefied_mean_bin_0.5/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.5

qiime tools export --input-path rarefied_mean_bin_0.6/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.6
qiime tools export --input-path rarefied_mean_bin_0.6/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.6
# 
# qiime tools export --input-path rarefied_mean_bin_0.7/feature_importance.qza --output-path ALL_rarefied_mean_bin/features_0.7
# qiime tools export --input-path rarefied_mean_bin_0.7/predictions.qza --output-path ALL_rarefied_mean_bin/predictions_0.7


  
