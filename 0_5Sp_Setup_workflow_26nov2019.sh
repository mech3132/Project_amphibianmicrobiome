#!bin/bash

# Paths to download from QIITA
## NOTE NEED TO UPDATE BECAUSE YOU CANNOT CURRENTLY DOWNLOAD THESE FILES BECAUSE I HAVEN'T MADE QIITA PUBLIC
featureSequence1=''
featureSequence2=''
featureSequence3=''
featureTable='' # NOT rarefied
mappingFile=''
betaUWU=''
alphaFaithPD=''
alphaObsOTU=''
tree=''
DOWNLOAD=FALSE
RUN_CLASSIFICATION=FALSE # For running sklearn classificaiton. Set to FALSE because it takes forever each time.

# Download into directory for all QIITA files
if $DOWNLOAD; then
    mkdir 0_QIITA_files
    cd 0_QIITA_files
    curl $featureTable -o ./0_QIITA_files/featureTable.qza
    curl $mappingFile -o ./0_QIITA_files/mapping_5sp.txt
    curl $tree -o ./0_QIITA_files/tree16_5Sp.tre
    cd ..
    
    mkdir 0_fasta_files
    cd 0_fasta_files
    curl -O $featureSequence1
    curl -O $featureSequence2
    curl -O $featureSequence3
    cd ..
fi

mkdir ./0_QIITA_files/intermediate_files

# First, fix mapping file because it has duplicate LinkerSequencePrimer columns
# Get column number of LinkerPrimerSequence
NUMBER=$(head -1 ./0_QIITA_files/mapping_5sp.txt | tr '\t' '\n' | cat -n | grep "LinkerPrimerSequence" | awk '{print $1}')
# Remove from mapping file
# The problem with awk is that it doesn't remove tabs. So, I'm doing a work-around where you replace it all with 'blank'
awk -v col=$NUMBER -F $'\t' 'BEGIN {OFS = FS} {$col="blank";print}' ./0_QIITA_files/mapping_5sp.txt > ./0_QIITA_files/intermediate_files/mapping_5sp_fixed.txt 
qiime tools inspect-metadata ./0_QIITA_files/intermediate_files/mapping_5sp_fixed.txt 


# Filter OTU table
qiime feature-table filter-samples \
--i-table ./0_QIITA_files/featureTable.qza \
--m-metadata-file ./0_QIITA_files/intermediate_files/mapping_5sp_fixed.txt \
--p-where "prepost_exposure IN ('Pre', 'Post')" \
--o-filtered-table ./0_QIITA_files/intermediate_files/featureTable_samplFilt.qza 



if $RUN_CLASSIFICATION ; then

	######## ASSIGNING TAXONOMY ############

	# Convert 3 fasta files into one
	cat ./0_fasta_files/79610_reference-hit.seqs.fa ./0_fasta_files/79613_reference-hit.seqs.fa ./0_fasta_files/79616_reference-hit.seqs.fa > ./0_fasta_files/allSeqs.fna
	# remove duplicates
	awk '!a[$0]++' ./0_fasta_files/allSeqs.fna > ./0_fasta_files/allSeqs_unique.fna
	rm ./0_fasta_files/allSeqs.fna

	# Count number of sequences "before"
	touch ./0_fasta_files/LOG
	echo "Number of sequences in allSeqs_unique.fna before filtering: " > ./0_fasta_files/LOG
	grep -c ">" ./0_fasta_files/allSeqs_unique.fna >> ./0_fasta_files/LOG
	echo "" >> ./0_fasta_files/LOG

	qiime tools import \
	--input-path ./0_fasta_files/allSeqs_unique.fna \
	--output-path ./0_fasta_files/allSeqs_unique.qza \
	--type FeatureData[Sequence]

	gzip ./0_fasta_files/*.fa
	rm ./0_fasta_files/*.fna

	# Now, filter fasta so that it only has OTUs belonging in otu table
	qiime feature-table filter-seqs \
	--i-data ./0_fasta_files/allSeqs_unique.qza \
	--i-table ./0_QIITA_files/intermediate_files/featureTable_samplFilt.qza \
	--o-filtered-data ./0_fasta_files/filtered_repset.qza

	# Count number of sequences "after"
	qiime tools export \
	--input-path ./0_fasta_files/filtered_repset.qza \
	--output-path ./0_fasta_files/temp_filtered_repset

	echo "Number of sequences in filtered_repset.qza before filtering: " >> ./0_fasta_files/LOG
	grep -c ">" ./0_fasta_files/temp_filtered_repset/dna-sequences.fasta >> ./0_fasta_files/LOG
	echo "" >> ./0_fasta_files/LOG


	# Clean up
	rm ./0_fasta_files/allSeqs_unique.qza
	rm -rf ./0_fasta_files/temp_filtered_repset

	# Just so you know, the "before" filtering had 16008 seqs, and the after filtering had 15000 seqs

	# Going to try matching to 2 different databases. See which one has more hits.

	## I tried doing gg and silva; silva had more hits so I'm going to use that here
    ## Download greengenes taxonomy 
    curl -O "https://s3-us-west-2.amazonaws.com/qiime2-data/2019.10/common/gg-13-8-99-515-806-nb-classifier.qza"

    # Add taxonomy to featureTable
    qiime feature-classifier classify-sklearn \
     --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
     --i-reads ./0_fasta_files/filtered_repset.qza \
     --o-classification ./0_QIITA_files/taxonomy_gg.qza

    # Download silva taxonomy 
#     curl -O "https://s3-us-west-2.amazonaws.com/qiime2-data/2019.10/common/silva-132-99-515-806-nb-classifier.qza"

    # Add taxonomy to featureTable
#     qiime feature-classifier classify-sklearn \
#       --i-classifier silva-132-99-515-806-nb-classifier.qza \
#       --i-reads ./0_fasta_files/filtered_repset.qza \
#       --o-classification ./0_QIITA_files/taxonomy.qza
    # STOPPED RUNNING HERE TEMP 27nov2019
    # remove database to save space
#     rm silva-132-99-515-806-nb-classifier.qza
    rm gg-13-8-99-515-806-nb-classifier.qza

    # Rename taxa to either b gg or silva
    mv ./0_QIITA_files/taxonomy_gg.qza ./0_QIITA_files/taxonomy.qza
    
    

fi

# I did run both Greengenes and SILVA-- although SILVA performed better, I found that GG had more hits for inhibitory strains so I use that instead. Here were the number of hits using both:
# 
# SILVA:
# Kingdom: 15948
# Phylum: 15388
# Class: 15106
# Order: 14140
# Family: 12898
# Genera: 8919\
# 
# GG:
# Kingdom: 15982
# Phylum: 15266
# Class: 14496
# Order: 13553
# Family: 12109
# Genera: 8719\



####### ADD TAXONOMY TO TABLE 
# Export
qiime tools export \
--input-path ./0_QIITA_files/taxonomy.qza \
--output-path ./0_QIITA_files/intermediate_files/exported-taxonomy
# Change header names
sed '1d' ./0_QIITA_files/intermediate_files/exported-taxonomy/taxonomy.tsv > ./0_QIITA_files/intermediate_files/exported-taxonomy/taxonomy_temp.tsv
echo -e "#OTUID\ttaxonomy\tconfidence" | cat - ./0_QIITA_files/intermediate_files/exported-taxonomy/taxonomy_temp.tsv > ./0_QIITA_files/intermediate_files/exported_taxonomy_wheaders.tsv
rm -rf ./0_QIITA_files/intermediate_files/exported-taxonomy

# Export biom so we can load into R
qiime tools export --input-path ./0_QIITA_files/intermediate_files/featureTable_samplFilt.qza --output-path ./0_QIITA_files/intermediate_files/featureTable_count

# Add taxonomy to feature Table
biom add-metadata -i ./0_QIITA_files/intermediate_files/featureTable_count/feature-table.biom \
-o ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.biom \
--observation-metadata-fp ./0_QIITA_files/intermediate_files/exported_taxonomy_wheaders.tsv \
--sc-separated taxonomy \
--observation-header OTUID,taxonomy,confidence

# Check how many of these taxa are "unidentified" or "unknown"
biom convert -i ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.biom \
--to-tsv \
--header-key taxonomy \
-o ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt 

touch ./0_QIITA_files/LOG
echo "Number of Unassigned:" > ./0_QIITA_files/LOG
grep "Unassigned" -c ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt >> ./0_QIITA_files/LOG
echo "Number of Chloroplasts:" >> ./0_QIITA_files/LOG
grep "Chloroplast" -c ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt  >> ./0_QIITA_files/LOG
echo "Number of Mitochondria:" >> ./0_QIITA_files/LOG
grep "Mitochondria" -c ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt >> ./0_QIITA_files/LOG
echo "" >> ./0_QIITA_files/LOG

#Now, let's filter some sequences to remove mito, chloro, eukaryotes, etc

touch ./0_QIITA_files/intermediate_files/features_to_remove.txt 
echo "Unassigned" > ./0_QIITA_files/intermediate_files/features_to_remove.txt
echo "Chloroplast" >> ./0_QIITA_files/intermediate_files/features_to_remove.txt
echo "Mitochondria" >> ./0_QIITA_files/intermediate_files/features_to_remove.txt


# Run this, and then MANUALLY curate "contaminants"
RScript R_code/decontam_workflow_26nov2019.R

# Streptococcus
echo "TACGTAGGTCCCGAGCGTTGTCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAGTCTGAAGTTAAAGGCTGTGGCTTAAC" >> ./0_QIITA_files/intermediate_files/features_to_remove.txt

# Remove these
grep -v -f ./0_QIITA_files/intermediate_files/features_to_remove.txt ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa.txt >> \
 ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa_filtered.txt

# Now, move all files that I want to use for downstream analyses into a folder called "DOWNSTREAM"
mkdir ./0_QIITA_files/downstream
mv ./0_QIITA_files/intermediate_files/featureTable_count_wtaxa_filtered.txt ./0_QIITA_files/downstream
mv ./0_QIITA_files/intermediate_files/mapping_5sp_fixed.txt ./0_QIITA_files/downstream

#### Finally, convert to biom and qza for downstream analysis
# Convert and import OTU table
biom convert -i ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.txt \
--to-hdf5 \
--header-key taxonomy \
-o ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.biom

qiime tools import \
--input-path ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.biom \
--output-path ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.qza \
--type 'FeatureTable[Frequency]'


# Rarefy the OTU table
qiime feature-table rarefy \
--i-table ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.qza \
--p-sampling-depth 10000 \
--o-rarefied-table ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.qza

# Finally, convert to text so we can use it in downstream
qiime tools export \
--input-path ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.qza \
--output-path ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered

mv ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered/feature-table.biom ./0_QIITA_files/downstream/temp_notaxa.biom
rm -rf ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered

# We need to re-add metadata because qza doesn't have metadata attached to it.
biom add-metadata -i ./0_QIITA_files/downstream/temp_notaxa.biom \
-o ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.biom \
--observation-metadata-fp ./0_QIITA_files/intermediate_files/exported_taxonomy_wheaders.tsv \
--sc-separated taxonomy \
--observation-header OTUID,taxonomy,confidence

rm ./0_QIITA_files/downstream/temp_notaxa.biom

biom convert -i ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.biom \
--to-tsv \
--header-key taxonomy \
-o ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.txt

# Lastly, record feature table information
echo "" >> ./0_QIITA_files/LOG
echo "Rarefied (10k) table stats:" >> ./0_QIITA_files/LOG
biom summarize-table -i ./0_QIITA_files/downstream/featureTable_r10k_wtaxa_filtered.biom | head -n 4 >> ./0_QIITA_files/LOG
echo "" >> ./0_QIITA_files/LOG

echo "" >> ./0_QIITA_files/LOG
echo "NON-rarefied table stats:" >> ./0_QIITA_files/LOG
biom summarize-table -i ./0_QIITA_files/downstream/featureTable_count_wtaxa_filtered.biom >> ./0_QIITA_files/LOG
echo "" >> ./0_QIITA_files/LOG
