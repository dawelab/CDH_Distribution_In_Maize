## 1. Process data
1. Select only reads greater than 64bp. I obtained this data which has already been demultiplexed and had qualified reads selected. See Issues for more details. This resulted in reads that were as small as 22bp. TASSEL cannot take anything shorter than 64bp. I had to remove them. 
2. Run barcode faker on the Romero-Navarro size filtered data to make the data compatible with TASSEL. TASSEL provides this script for demultiplexing, I modified it to set the barcode length to 13 in order to have enough permutations of ACTG in for the number of lines that I have.
3. Modify barcode faker file names to be compatible with TASSEL and my other reads.
4. This generates the key file needed for TASSEL for all data sources. DC1 = Dawe Lab Controls 1, DC2 = Dawe Lab Controls 2, SW = Swarts et al 2017, RN = Romero-Navarro 2018, RY = Romay et al 2013. 
5. This is the final key file used
6. This appends the K10L2 contig to the Mo17 genome. 
7. This makes a file with all the metadata for all the datasets. 1.7_Files includes all the files that the script loads in.
