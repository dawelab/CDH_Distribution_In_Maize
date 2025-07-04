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
>5. This is the output file of 4
>6. This uses a random forest model (ntrees=1000000) to assign Ab10 class. It requires that 65% of the decision trees report the same class and performs a PCA on the highest mean decreasing gini score variables (=2) to explore the ambigous classes.
>7. This writes out the Ab10 classes and generates a few tables describing what Maize type Ab10 and the types occur in.
>8. This plots all of the Ab10 and N10 controls along with all of the Ab10 called positive experimental samples



3. This determines if all samples that were identified as N10 above are actually N10 or K10L2
>1. This uses Kmeans clusters to determine if controls are N10 or K10L2. I am separating the heterozygous (1) and homozygous (2) controls because it improves the clustering.
>2. This uses Kmeans clusters to class all experimental samples called N10 in 6.2.2 into K10L2 and N10 100 times. It requires that all the control samples be classed correctly.
>3. This script merges all the individual files into a single output file. I ran it on the command line. 
>4. This is the final call file edited so that Ab10 = 1 and N10=0.
>5. This uses Kmeans clusters to determine if controls are N10 or K10L2. I am including Ab10-I and Ab10-III in the model because it improves clutering. I am separating the heterozygous and homozygous (2) controls because it improves the clustering.
>6. This uses Kmeans clusters to class all experimental samples called N10 in 6.2.2 into K10L2 and N10 100 times. I am including Ab10-I and Ab10-III in the model because it improves clutering. It requires that all the control samples be classed correctly.
>7. This is the final call file edited so that K10L2 = 1 and N10=0.
>8. This writes out the K10L2 calls in the group file.
>9. This plots all of the K10L2 and N10 controls along with all of the K10L2 called positive experimental samples

4. This determines if all samples have B Chromosomes 
>1. This creates files used to annotate the B Chromosome.
>2. This uses Kmeans clusters to determine if controls do or do not have B chromosomes. I am separating the high (1) and low (2) copy number controls because it improves the clustering. The low B chromosome copy number controls could not be subsampled into 3 groups.
>3.  This uses Kmeans clusters to class all experimental samples with the high copy number B chromosomes. It requires that all the control samples be classed correctly. This must happen first or the controls will not be reliably classed and the script will hang.
>4.  This joins the Kmeans clusters scripts and the quality control scripts.
>5.   This uses Kmeans clusters to class all experimental samples that were not called as B positive in 6.4.3 with the low copy number B chromosome controls.
>6.   This joins the Kmeans clusters scripts and the quality control scripts.
>7.   This plots all the B chromosome samples identified and writes out a file with all the B Chromosome Positive Samples.
>8.   This makes the final calls for B chromosome positive, ambigious, and absent and writes them into the Groups File

5.   This divides each bin across the selfish genetic element by the mean of all single copy core gene tag index values. This acts as a proxy for copy number for the B Chromosome validated by Abnormal Chromsoome 10.
>1. Ab10
>2. K10L2
>3. BChrom
>4. Plot all copy numbers together
>5. Merge all of the copy number data with the large meta deta file.

6. This is the groups file with the final calls for Ab10, K10L2, and the B chromosome including pseudo copy number
