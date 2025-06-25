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
