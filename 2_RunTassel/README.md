These scripts run the full [TASSEL pipeline](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/Tassel5GBSv2Pipeline/GBSSeqToTagDBPlugin) for the B73 Ab10 and B chromosome reference and the CI66 K10L2 assembly. The Ab10 and K10L2 haplotype assemblies are available through https://doi.org/10.1093/genetics/iyaf091. The B chromosome assembly is from https://doi.org/10.1073/pnas.2104254118. 

## 2. Run TASSEL
1. This runs everything for TASSEL. It ran in less than a day and took about 400GB of memory.
> 1. B73 Ab10 and B Chr. 
> 2. K10L2

2. This subsets the big Tag by Taxa file down to the Ab10 and B chromosome.
> 1. This identifies the Ab10 and B chromosome Tags l
> 2. This is the table used in the array job for 2.2.3
> 3. This extracts Ab10 and the B chromosome from each Tassel run.

3. This subsets the big Tag by Taxa file down to K10L2
> 1. This identifies K10L2 Tags for each run of TASSEL.
> 2. This is the table used in the array job for 2.3.3
> 3. This extracts K10L2 from each Tassel run.
