#This defines variables
DIR=""
DIRNAME=""
TAXA="Tassel_TagTaxaDist_AllData_v7_v_K10L2"

join $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.15.txt $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.14.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.1
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.1 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.13.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.2
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.2 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.12.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.3
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.3 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.11.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.4
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.4 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.10.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.5
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.5 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.9.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.6
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.6 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.8.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.7
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.7 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.7.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.8
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.8 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.6.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.9
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.9 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.5.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.10
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.10 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.4.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.11
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.11 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.3.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.12
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.12 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.2.txt > $DIR/$DIRNAME/Subset_TaxaFiles/temp.13
join $DIR/$DIRNAME/Subset_TaxaFiles/temp.13 $DIR/$DIRNAME/Subset_TaxaFiles/$TAXA.Sub1Perc.Missing.1.txt > $DIR/$DIRNAME/$TAXA.Sub1Perc.Missing.txt

rm $DIR/$DIRNAME/Subset_TaxaFiles/temp.*
