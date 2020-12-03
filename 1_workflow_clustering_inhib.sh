#!bin/bash

mkdir 1_woodhams_antifungal
cd 1_woodhams_antifungal

##### File paths ####
woodhams_seq="http://www.esapubs.org/archive/ecol/E096/059/Amphibian-skin_bacteria_16S_sequences.fna"
woodhams_meta="http://www.esapubs.org/archive/ecol/E096/059/Amphibian-skin_bacteria_metadata.txt"
featureTable='../0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.qza'
featureSeq='../0_fasta_files/filtered_repset.qza'

##### Data import and filtering
# First, download the amphibian database from woodhams
curl -O $woodhams_seq
curl -O $woodhams_meta
# BUT line 1138 had a random sequence with ":" in it, so I just deleted that one.
grep -B1 ":" Amphibian-skin_bacteria_16S_sequences.fna > omit.txt
grep -v -f omit.txt Amphibian-skin_bacteria_16S_sequences.fna > Amphibian-skin_bacteria_16S_sequences2.fna 
rm Amphibian-skin_bacteria_16S_sequences.fna
rm omit.txt

qiime tools import \
--input-path Amphibian-skin_bacteria_16S_sequences2.fna \
--output-path Amphibian-skin_bacteria_16S_sequences.qza \
--type 'FeatureData[Sequence]'

# Filter sequence data 
qiime feature-table filter-seqs \
--i-data $featureSeq \
--i-table $featureTable \
--o-filtered-data final_filtered_repset.qza

##### VSEARCH #####
# Finally, I ran vsearch in qiime2

qiime vsearch cluster-features-closed-reference \
--i-sequences final_filtered_repset.qza \
--i-table $featureTable \
--i-reference-sequences Amphibian-skin_bacteria_16S_sequences.qza \
--p-perc-identity 1 \
--output-dir hits_inhib_database

# Finally, I exported the files so I can use them easier
qiime tools export \
--input-path ./hits_inhib_database/clustered_sequences.qza \
--output-path ./hits_inhib_database/exported_clustered_sequences

# First, get list of seq
grep -A1 'inhibitory' ./hits_inhib_database/exported_clustered_sequences/dna-sequences.fasta \
| grep -v 'inhibitory' | grep -v '\-\-'> ./hits_inhib_database/exported_clustered_sequences/allSeqs.txt
# Then, get list of inhibitory
grep 'inhibitory' ./hits_inhib_database/exported_clustered_sequences/dna-sequences.fasta | sed 's/>//g' > ./hits_inhib_database/exported_clustered_sequences/otu_names.txt

paste ./hits_inhib_database/exported_clustered_sequences/otu_names.txt ./hits_inhib_database/exported_clustered_sequences/allSeqs.txt > mapping_inhib_otus.txt

echo "otu_id\tSequence" | cat - mapping_inhib_otus.txt > mapping_inhib_otus_final.txt
