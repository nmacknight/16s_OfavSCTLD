# Qiime2 Pipeline :microbe:
**Project Context:** Approximately 2,500 coral fragments of the reef-building coral, Orbicella faveolata comprising 179 unique genotypes were exposed to Stony Coral Tissue Loss Disease (SCTLD) in 2021 at Mote Marine Laboratory for a Florida DEP Funded Project.

Here, we process ~1650 of those samples through the Qiime2 pipeline for 16s microbial ecology analysis.
Leveraging 192 unique barcodes, we performed the library prep in house and sequenced the samples across 9 sequencing runs (192 samples/run). 

> Sequenced the 515F to 806R V4 region of the 16s rRNA on an Illumina MiSeq v3 PE300 with a target read count of 20 Million/run or ~100,000/sample.

**Note:** Steps 1-3 below was applied to each sequencing run. Rather than combine all 9 runs into one script, I provided simply the first run ("D1D2") to demonstrate how each run was processed. The analysis steps from 1-3 are the same for each run, and then each run is merged together at step 4. 

## Sequencing run - "D1D2"

```
cd Desktop/Ofav\ SCTLD/Raw\ Data/D1D2-230117_M02476_0565_000000000-KRVDW/Alignment_1/20230120_001152

mkdir ./DemultiplexedSeqs
mv ./Fastq/FastqSummaryF1L1.txt ./../
mv ./Fastq/Undetermined* ./../
```
So I removed the FastqSummary and Undetermined (which are reads whose sample they belong to are undetermined). Because they were stored in the Fastq folder but shouldnt actually be in there because they are not our samples Fastq files. 

# 1 Import Data

With this style of import, the samples are ALREADY Demultiplexed.

Importing Data - **Casava 1.8 paired-end demultiplexed fastq**

**Note:** You will receive files in different formats from your sequencer. It may not be exactly like mine. For example in my White Plague Qiime2 folder you can see they are different. Such as if your reads are single end v paired end reads, multiplexed or not, or demultiplexed or still multiplexed. 
Mine are paired end demultiplexed reads that were received in a format that reflects the Casava 1.8 paired-end demultiplexed fastq style. 
To Read up on which is best for your project checkout Qiime2's [Importing Data](https://docs.qiime2.org/2023.2/tutorials/importing/) page.

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Fastq \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza
```
>Run Time: 20 min

```
mv demux-paired-end.qza ./DemultiplexedSeqs/
```

Summary of the demultiplexing results

```
qiime demux summarize \
  --i-data ./DemultiplexedSeqs/demux-paired-end.qza \
  --o-visualization ./DemultiplexedSeqs/demux-paired-end.qzv
```
> RUN TIME: 12 Min.

View The Interactive Plot. Truncate when black bars begin to consistently go below a Quality Score of 30. 

```
qiime tools view ./DemultiplexedSeqs/demux-paired-end.qzv
```
> RUN TIME: <1 min.

We will need to remove the adapters before the reads are truncated. To do this, we use Cutadapt in step 2.


# 2. Cutadapt

**Removes Primer from Reads**

This trim-paired form of cutadapt will remove any nucelotides that are before and including the provided sequence on the front and reverse read of the sequence. This removes non-biological sequence data. 

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences ./DemultiplexedSeqs/demux-paired-end.qza \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r GGACTACNVGGGTWTCTAAT \
  --o-trimmed-sequences ./DemultiplexedSeqs/demux-paired-end-trimmed.qza
```

> Run time: (forgot to track this one, dont be mad, its good time to get a snack.)

Summary of the demultiplexing results

```
qiime demux summarize \
  --i-data ./DemultiplexedSeqs/demux-paired-end-trimmed.qza  \
  --o-visualization ./DemultiplexedSeqs/demux-paired-end-trimmed.qzv
```

>RUN TIME: 12 Min.

View The Interactive Plot. Truncate when black bars begin to consistently go below a Quality Score of 30.

```
qiime tools view ./DemultiplexedSeqs/demux-paired-end-trimmed.qzv
```



# 3. DADA2

**Chops the read length where the quality begins to become poor.**
**Note:** Update the number after "--p-trunc-len-f" and "--p-trunc-len-r" to reflect YOUR cutoff point as determined after looking at the figures produced by the command above.

```
mkdir ./DADA2
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./DemultiplexedSeqs/demux-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trunc-len-f 280 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 204 \
  --o-representative-sequences ./DADA2/rep-seqs-dada2.qza \
  --o-table ./DADA2/table-dada2.qza \
  --o-denoising-stats ./DADA2/stats-dada2.qza
```
> RUN TIME: 40 minutes for (463mb demux.qza)), or 14 hours (11gb demux.qza file). File size depends on read depth and number of samples you are working with. 

```
qiime metadata tabulate \
  --m-input-file ./DADA2/stats-dada2.qza \
  --o-visualization ./DADA2/stats-dada2.qzv
qiime feature-table summarize \
  --i-table ./DADA2/table-dada2.qza \
  --o-visualization ./DADA2/table-dada2.qzv \
  --m-sample-metadata-file D1D2-metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data ./DADA2/rep-seqs-dada2.qza \
  --o-visualization ./DADA2/rep-seqs-dada2.qzv
```
> Run Time < 1 min.

### Now Repeat these steps above for each sequencing run before merging below 

#From her we are merging the data 
#For now we have the following datasets analyzed in Qiime2

**Note:** Below are the directories of where the results are stored of running the above code for each specific sequencing run. 

/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D1D2-230117_M02476_0565_000000000-KRVDW/Alignment_1/20230120_001152/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D3D4-230120_M02476_0566_000000000-KT3CN/Alignment_1/20230123_004935/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D5D6-230123_M02476_0567_000000000-KT3F4/Alignment_1/20230125_225417/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D7D8-230127_M02476_0568_000000000-KRVYJ/Alignment_1/20230129_230519/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D9D18-230131_M02476_0569_000000000-KT4LR/Alignment_1/20230203_005408/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D11D10-230206_M02476_0570_000000000-KRNHM/Alignment_1/20230209_022123/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D12D13-230213_M02476_0571_000000000-KRVDL/Alignment_1/20230216_020440/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D14D15-230217_M02476_0572_000000000-KT4LV/Alignment_1/20230219_235336/DADA2/
/Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/D16D17-230221_M02476_0573_000000000-KT3DP/Alignment_1/20230224_005823/DADA2/


# 3. Merging Sequencing Runs

Making a folder for all this work to go into:

```
cd /Users/nicholas.macknight/Desktop/Ofav SCTLD/Raw Data/
mkdir Merged
cd Merged
```

```
qiime feature-table merge \
  --i-tables ../D1D2-230117_M02476_0565_000000000-KRVDW/Alignment_1/20230120_001152/DADA2/table-dada2.qza \
  --i-tables ../D3D4-230120_M02476_0566_000000000-KT3CN/Alignment_1/20230123_004935/DADA2/table-dada2.qza \
  --i-tables ../D5D6-230123_M02476_0567_000000000-KT3F4/Alignment_1/20230125_225417/DADA2/table-dada2.qza \
  --i-tables ../D7D8-230127_M02476_0568_000000000-KRVYJ/Alignment_1/20230129_230519/DADA2/table-dada2.qza \
  --i-tables ../D9D18-230131_M02476_0569_000000000-KT4LR/Alignment_1/20230203_005408/DADA2/table-dada2.qza \
  --i-tables ../D11D10-230206_M02476_0570_000000000-KRNHM/Alignment_1/20230209_022123/DADA2/table-dada2.qza \
  --i-tables ../D12D13-230213_M02476_0571_000000000-KRVDL/Alignment_1/20230216_020440/DADA2/table-dada2.qza \
  --i-tables ../D14D15-230217_M02476_0572_000000000-KT4LV/Alignment_1/20230219_235336/DADA2/table-dada2.qza \
  --i-tables ../D16D17-230221_M02476_0573_000000000-KT3DP/Alignment_1/20230224_005823/DADA2/table-dada2.qza \
  --o-merged-table merged-table.qza
```
>Run time: < 1 min.

```
qiime feature-table merge-seqs \
  --i-data ../D1D2-230117_M02476_0565_000000000-KRVDW/Alignment_1/20230120_001152/DADA2/rep-seqs-dada2.qza \
  --i-data ../D3D4-230120_M02476_0566_000000000-KT3CN/Alignment_1/20230123_004935/DADA2/rep-seqs-dada2.qza \
  --i-data ../D5D6-230123_M02476_0567_000000000-KT3F4/Alignment_1/20230125_225417/DADA2/rep-seqs-dada2.qza \
  --i-data ../D7D8-230127_M02476_0568_000000000-KRVYJ/Alignment_1/20230129_230519/DADA2/rep-seqs-dada2.qza \
  --i-data ../D9D18-230131_M02476_0569_000000000-KT4LR/Alignment_1/20230203_005408/DADA2/rep-seqs-dada2.qza \
  --i-data ../D11D10-230206_M02476_0570_000000000-KRNHM/Alignment_1/20230209_022123/DADA2/rep-seqs-dada2.qza \
  --i-data ../D12D13-230213_M02476_0571_000000000-KRVDL/Alignment_1/20230216_020440/DADA2/rep-seqs-dada2.qza \
  --i-data ../D14D15-230217_M02476_0572_000000000-KT4LV/Alignment_1/20230219_235336/DADA2/rep-seqs-dada2.qza \
  --i-data ../D16D17-230221_M02476_0573_000000000-KT3DP/Alignment_1/20230224_005823/DADA2/rep-seqs-dada2.qza \
  --o-merged-data merged-rep-seqs.qza
```
>Run time: < 1 min.


Summary of merged artifact. You will need to concatenate the other sequencing run metadata files to create "sample-metadata.tsv" in excel/Numbers. I did it in Numbers (the only time ive worked in this program over excel) because Numbers can export by tsv. 


```
qiime feature-table summarize \
  --i-table merged-table.qza \
  --o-visualization merged-table.qzv \
  --m-sample-metadata-file sample-metadata.tsv
```
>Run time: < 1 min.

View
```
qiime tools view merged-table.qzv
```

```
qiime feature-table tabulate-seqs \
	--i-data merged-rep-seqs.qza \
	--o-visualization merged-rep-seqs.qzv
```
>Run Time <1 min.


# 4. TAXONOMIC CLASSIFICATION

```
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515-806-nb-classifier.qza \
  --i-reads merged-rep-seqs.qza \
  --o-classification taxonomy.qza
```
> Run time ~17 hours (merged-rep-seqs.qza was 12.7mb)

```
qiime feature-table filter-features \
  --i-table merged-table.qza \
  --m-metadata-file taxonomy.qza \
  --o-filtered-table id-filtered-table.qza
```
>Run time: < 1 min.

```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```
>Run Time: < 1 min.

```
qiime tools view taxonomy.qzv
```

```
qiime taxa barplot \
  --i-table merged-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```
>Run Time: ~ 2 min.

```
qiime tools view taxa-bar-plots.qzv
```

**Filtering out Chloroplast and Mitochondria From Dataset**
There is a lot of chloroplast from Osteobium sp. perhaps. Lets remove them.

```
qiime taxa filter-table \
  --i-table merged-table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --p-include d__Bacteria,d__Archaea \
  --o-filtered-table table-BacArc.qza
```
>Run Time: < 1 min.

```
qiime taxa filter-seqs \
  --i-sequences merged-rep-seqs.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --p-include d__Bacteria,d__Archaea \
  --o-filtered-sequences rep_seqs_BacArc.qza
```
>Run Time: < 1 min.

```
qiime feature-table summarize \
  --i-table table-BacArc.qza \
  --o-visualization table-BacArc.qzv \
  --m-sample-metadata-file sample-metadata.tsv
```
>Run Time: < 1 min.

```
qiime feature-table tabulate-seqs \
--i-data rep_seqs_BacArc.qza \
--o-visualization rep_seqs_BacArc.qzv
```
>Run Time: < 1 min.

```
qiime taxa barplot \
  --i-table table-BacArc.qza  \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv\
  --o-visualization taxa-bar-BacArc.qzv
```
>Run Time: ~ 2 min.

```
qiime tools view taxa-bar-BacArc.qzv
```


# 5. GENERATE A TREE FOR PHYLOGENETIC DIVERSITY ANALYSES

```
mkdir PhylogeneticTree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep_seqs_BacArc.qza \
  --o-alignment ./PhylogeneticTree/aligned-rep-seqs.qza \
  --o-masked-alignment ./PhylogeneticTree/masked-aligned-rep-seqs.qza \
  --o-tree ./PhylogeneticTree/unrooted-tree.qza \
  --o-rooted-tree ./PhylogeneticTree/rooted-tree.qza
```
>RunTime: ~2.5 days!


# 6. Export

```
qiime tools export \
	--input-path table-BacArc.qza \
	--output-path table
  
qiime tools export \
	--input-path taxonomy.qza \
	--output-path taxonomy
```
>Run Time: 1 min.

```
qiime tools export \
	--input-path ./PhylogeneticTree/rooted-tree.qza \
	--output-path rooted-tree
```
>Run Time: < 1 min.

Congrats, ya did it! :clinking_glasses:
