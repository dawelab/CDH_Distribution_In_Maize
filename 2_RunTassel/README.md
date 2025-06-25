## 2. Run TASSEL
1. This runs everything for TASSEL. It ran in less than a day and took about 400GB of memory.
> 1. This runs it using the corrected version of the HiFi B73-Ab10 version from Mingyu.
> 2. These scripts are aligning to the K10L2 assembly. I'm testing if the unscaffolded contigs or the the K10L2 contig and the Mo17 genome work better.

2. This subsets the big Tag by Taxa file down to the Ab10 and B chromosome.
> 1. This identifies the Ab10 and B chromosome Tags l
> 2. This is the table used in the array job for 2.2.3
> 3. This extracts Ab10 and the B chromosome from each Tassel run.

3. This subsets the big Tag by Taxa file down to K10L2
> 1. This identifies K10L2 Tags for each run of TASSEL.
> 2. This is the table used in the array job for 2.3.3
> 3. This extracts K10L2 from each Tassel run.
