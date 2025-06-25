## 3. Normalize Tag Read Depth
### Ab10
1. This subsets the large Tag by Taxa File from TASSEL down to 1% 15 times.
2. This calculates the sum of all mapped tags in each randomly selected 1% of the genome. 
3. This merges all of the sum files for each of the 15 1% sub-samples of the genome. 
4. This normalizes the coverage on both Ab10 and the B chromosome by the average sum of 1% of tags per sample for all three runs.

### K10L2
1. This subsets the large Tag by Taxa File from TASSEL down to 1% 15 times.
2. This calculates the sum of all mapped tags in each randomly selected 1% of the genome. It uses the same script as the Ab10 3.2, but applied to the CI66-K10L2 files.
3. This merges all of the sum files for each of the 15 1% sub-samples of the genome.
4. This normalizes the coverage on both Ab10 and the B chromosome by the average sum of 1% of tags per sample for all three runs. It uses the same script as Ab10 3.4, but applied to the CI66-K10L2 files. 
