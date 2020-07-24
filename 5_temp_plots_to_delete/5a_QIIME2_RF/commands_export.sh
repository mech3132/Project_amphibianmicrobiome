

# Note to self: export predictions and feature_importance to plot in R manually.

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

