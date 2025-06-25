## 3. Normalize Tag Count
### Ab10
1. This subsets the large Tag by Taxa File down to 1% 15 times for each TASSEL run.
2. This calculates the sum of all mapped tags in the randomly selected 1% of the genome for each run of TASSEL. 
3. This merges all of the Sum files for each of the 15 1% sub-samples for each run of TASSEL. 
4. This normalizes the coverage on Ab10 and the B chromosome by the average sum of 1% of tags per sample for all three runs.

### K10L2
1. This subsets the large Tag by Taxa File down to 1% 15 times for each TASSEL run.
2. This calculates the sum of all mapped tags in the randomly selected 1% of the genome for each run of TASSEL. It uses the same R script in 3.2, but just feeds it differnt files.
3. This merges all of the Sum files for each of the 15 1% sub-samples for each run of TASSEL.
4. This normalizes the coverage on K10L2 by the average sum of 1% of tags per sample for all three runs. It uses the same script from 3.4, but just feeds it different files. 

## 4. Plot the Control Data
1. This scripts plots all the controls on the Ab10 haplotype.
2. This sums over Ab10 regions of known importance. It both does and does not include K10L2.
3. This plots all the controls for the B chromosome
4. This sums over the B chromosome.
5. This plots all the controls for K10L2
>1. This uses the K10L2 contigs
>2. This uses K10L2 appended to the chopped Mo17 genome (only one copy of shared region. 

## 5. Check Data Quality
1. This checks that the filtering by mapping quality is treating all my datasets roughly the same for all SGEs.
>1.  Ab10
>2. B chromosome
>3. K10L2


2. These scriripts filter out samples with too much missing data. 
>1.  This script uses the random 1 percent subsets generated in 3.1 to determine percent missing data for each sample in the B73-Ab10 hifi genome. By happenstance it also ran the other two TASSEL runs. Three samples from the Mo17 genomes failed for an unknown reason. All of the B73-Ab10 samples completed as expected.
>2.  This joins all the files only for the B73 HiFi Ab10 run of TASSEL
>3.  This plots out the missing data and filters it. The Romay and Romero-Navarro data sets have distributions that are a bit more different from my controls than I would like.
>4.  This script uses the random 1 percent subsets generated in 3a.1 to determine percent missing data for each sample in the K10L2 contigs.
>5.   This joins all the files for the K10L2 Contigs
>6.   This plots out the missing data and filters it for K10L2.

3. This pulls all the tags overlapping all of the single copy core genes for B73-Ab10.
>1. This extracts the tag sequence of all the single copy core genes and then subsets this list into 20 files.
>2. This finds the index of each tag in the large tag by taxa file, and pulls that index from the large tag by taxa file for each of the 20 subsets.
>3. This takes those indexes and extracts the appropriate line from the large tag by taxa file for each of the 20 subsets
>4. This brings together the single copy core gene tag by taxa entries from all 20 subsets.
>5. This generates the bin file for the single copy core genes. It produces bins that have 1MB of single copy core gene sequence wihtout regard to the actual physical distance covered. This is acceptable because the tags have already been filtered to only those overlapign single copy core genes.
>6. This adds bin values to th e 20 subset files generated in 5.3.3. Performing this on the merged file made in 5.3.4 takes too long.
>7. This loads in all the files generated in 5.3.6 and generates the tag index

4. This annotates the HiFi Ab10 genome
   
## 6. Classify All Data
1. These scripts filter, prepare, and generate the tag index so that clustering can be done.
>1. This script filters the data and merges the Romero Navarro technical replicates into a single column per biological sample for Ab10 and the B chromosome.
>2. This script takes the processed data and generates the Ab10 tag Index from it.
>3.  This script filters the data and merges the Romero Navarro technical replicates into a single column per biological sample for K10L2.
>4.  This script takes the processed data and generates the K10L2 tag Index from it.
>5.  This script filters the data and merges the Romero Navarro technical replicates into a single column per biological sample for B Chrom.
>6.  This script takes the processed data and generates the B chromosome tag Index from it.

2. This derermines if all samples are Ab10 or N10 and then identifies Ab10 type
>1. This uses Kmeans clusters to determine if controls are Ab10 or N10. I am separating the heterozygous (1) and homozygous (2) controls because it improves the clustering. 
>2. This uses Kmeans clusters to class all experimental samples into Ab10 and N10 100 times. It requires that all the control samples be classed correctly. 
>3. This script merges all the individual files into a single output file. I ran it on the command line.
>4. This is the final call file edited so that Ab10 = 1 and N10=0.
>6. This uses a random forest model (ntrees=1000000) to assign Ab10 class. It requires that 50% of the decision trees report the same class.
>7. This uses a random forest model (ntrees=1000000) to assign Ab10 class. It requires that 65% of the decision trees report the same class and performs a PCA on the highest mean decreasing gini score variables (=2) to explore the ambigous classes.
>8. This writes out the Ab10 classes and generates a few tables describing what Maize type Ab10 and the types occur in. 

3. This determines if all samples that were identified as N10 above are actually N10 or K10L2
>1. This uses Kmeans clusters to determine if controls are N10 or K10L2. I am separating the heterozygous (1) and homozygous (2) controls because it improves the clustering.
>2. This uses Kmeans clusters to class all experimental samples called N10 in 7.2.2 into K10L2 and N10 100 times. It requires that all the control samples be classed correctly.
>3. This script merges all the individual files into a single output file. I ran it on the command line. 218 K10L2 Positive, 38 Ambiguous. Cutoffs 94 for K10L2 positive, 13 for N10, everything inbetween ambiguous determined by control clustering. 
>4. This is the final call file edited so that Ab10 = 1 and N10=0.
>5. This uses Kmeans clusters to determine if controls are N10 or K10L2. I am including Ab10-I and Ab10-III in the model because it improves clutering. Ab10-II muddles it because it seems to be intermediate. I am only clustering on region 7:10 as they are the most consistenly covered in K10L2 and Ab10. I am separating the heterozygous and homozygous (2) controls because it improves the clustering.
>6. This uses Kmeans clusters to class all experimental samples called N10 in 7.2.2 into K10L2 and N10 100 times. I am including Ab10-I and Ab10-III in the model because it improves clutering. It requires that all the control samples be classed correctly.
>7. This is the final call file edited so that K10L2 = 1 and N10=0.
>8. This writes out the K10L2 calls in the group file.
>9. This plots all of the K10L2 and N10 controls along with all of the K10L2 called positive experimental samples

4. This determines if all samples have B Chromosomes 
>1. This creates files used to annotate the B Chromosome.
>2. This uses Kmeans clusters to determine if controls do or do not have B chromosomes. I am separating the high (1) and low (2) copy number controls because it improves the clustering. The low B chromosome copy number controls could not be subsampled into 3 groups because I didn't have enough low (only 4).
>3.  This uses Kmeans clusters to class all experimental samples with the high copy number B chromosomes. It requires that all the control samples be classed correctly. This must happen first or the controls will not be reliably classed and the script will hang. I tried doing the low copy Bs first and doing mixed controls, both hung.
>4.  This joins the Kmeans clusters scripts and the quality control scripts.
>5.   This uses Kmeans clusters to class all experimental samples that were not called as B positive in 7.4.3 with the low copy number B chromosomes.
>6.   This joins the Kmeans clusters scripts and the quality control scripts.
>7.   This plots all the B chromosome samples identified and writes out a file with all the B Chromosome Positive Samples.
>8.   This makes the final calls for B chromosome positive, ambigious, and absent and writes them into the Groups File
>9.   This divides each bin across the selfish genetic element by the mean of all single copy core gene tag index values. This acts as a proxy for copy number for the B Chromosome validated by Abnormal Chromsoome 10.
>>1. Ab10
>>2. K10L2
>>3. BChrom

5. This is the groups file with the final calls for Ab10, K10L2, and the B chromosome including pseudo copy number

## 7. Map 
1. This maps all the samples
