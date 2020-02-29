#!bin/bash

# First, set up file paths

featureTable='./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.qza'
metadata='./QIITA_files/downstream/mapping_5sp_fixed.txt'
tree16='./0_QIITA_files/tree16_5Sp.tre'

# Make directory
mkdir 2_diversity_outputs

# Import rooted tree
qiime tools import \
	--input-path $tree16 \
	--output-path ./2_diversity_outputs/rooted_tree.qza \
	--type 'Phylogeny[Rooted]'

#### ALPHA DIVERSITY #####
# List of alpha diversity metrics to use
touch ./2_diversity_outputs/alpha_metrics.txt
echo "observed_otus" > ./2_diversity_outputs/alpha_metrics.txt
echo "chao1" >> ./2_diversity_outputs/alpha_metrics.txt
echo "shannon" >> ./2_diversity_outputs/alpha_metrics.txt
touch ./2_diversity_outputs/alphaPD_metrics.txt
echo "faith_pd" > ./2_diversity_outputs/alphaPD_metrics.txt

while read alpha; do
	qiime diversity alpha \
	--i-table $featureTable \
	--p-metric ${alpha} --o-alpha-diversity ./2_diversity_outputs/${alpha}.qza 
	
	qiime tools export \
	--input-path ./2_diversity_outputs/${alpha}.qza \
	--output-path ./2_diversity_outputs/exported_${alpha}
done < ./2_diversity_outputs/alpha_metrics.txt

while read alpha; do
	qiime diversity alpha-phylogenetic \
	--i-table $featureTable \
	--i-phylogeny ./2_diversity_outputs/rooted_tree.qza \
	--p-metric ${alpha} --o-alpha-diversity ./2_diversity_outputs/${alpha}.qza 
	
	qiime tools export \
	--input-path ./2_diversity_outputs/${alpha}.qza \
	--output-path ./2_diversity_outputs/exported_${alpha}
done < ./2_diversity_outputs/alphaPD_metrics.txt


#### BETA DIVERSITY ####
# List of alpha diversity metrics to use
touch ./2_diversity_outputs/beta_metrics.txt
echo "braycurtis" > ./2_diversity_outputs/beta_metrics.txt
touch ./2_diversity_outputs/betaPD_metrics.txt
echo "unweighted_unifrac" > ./2_diversity_outputs/betaPD_metrics.txt
echo "weighted_unifrac" >> ./2_diversity_outputs/betaPD_metrics.txt

while read beta; do
	qiime diversity beta \
	--i-table $featureTable \
	--p-metric ${beta} --o-distance-matrix ./2_diversity_outputs/${beta}.qza 
	
	qiime tools export \
	--input-path ./2_diversity_outputs/${beta}.qza \
	--output-path ./2_diversity_outputs/exported_${beta}
done < ./2_diversity_outputs/beta_metrics.txt

while read beta; do
	qiime diversity beta-phylogenetic \
	--i-table $featureTable \
	--i-phylogeny ./2_diversity_outputs/rooted_tree.qza \
	--p-metric ${beta} --o-distance-matrix ./2_diversity_outputs/${beta}.qza 
	
	qiime tools export \
	--input-path ./2_diversity_outputs/${beta}.qza \
	--output-path ./2_diversity_outputs/exported_${beta}
done < ./2_diversity_outputs/betaPD_metrics.txt



# CITATIONS
# *Chao, A. (1984). “Non-parametric estimation of the number of classes in a population”.
# Faith. D.P. (1992). “Conservation evaluation and phylogenetic diversity”. Biological Conservation. (61) 1-10.
# DeSantis, T.Z., Hugenholtz, P., Larsen, N., Rojas, M., Brodie, E.L., Keller, K. Huber, T., Davis, D., Hu, P., Andersen, G.L. (2006). “Greengenes, a Chimera-Checked 16S rRNA Gene Database and Workbench Compatible with ARB”. Applied and Environmental Microbiology (72): 5069–5072.
# Shannon, C.E. and Weaver, W. (1949). “The mathematical theory of communication”. University of Illonois Press, Champaign, Illonois.
# Sorenson, T. (1948) “A method of establishing groups of equal amplitude in plant sociology based on similarity of species content.” Kongelige Danske Videnskabernes Selskab 5.1-34: 4-7.
# Lozupone, C. and Knight, R. (2005). “UniFrac: a new phylogenetic method for comparing microbial communities.” Applied and environmental microbiology 71 (12): 8228-8235.

