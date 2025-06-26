## 5. Check Data Quality
1. This checks that the filtering by mapping quality is treating all my datasets roughly the same for all SGEs. Plots derived from these are available in 5.1_Plots
>1.  Ab10
>2. B chromosome
>3. K10L2


2. These scriripts filter out samples with too much missing data. 
>1.  This script uses the random 1 percent subsets generated in 3.1 to determine percent missing data for each sample. These scripts use the same R script on different data sets for the B73 Ab10 B chr reference and the K10L2 reference.
>2.  This joins all the output files ofrom 2.1 for B73 Ab10  BChrom.
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
   
