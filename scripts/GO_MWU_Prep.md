## Preparing GO list and Results Files for [GO_MWU](https://github.com/z0on/GO_MWU) Analysis


GO_MWU requires a file with one column of gene names, and one column of the -log signed p-values (See [this](https://github.com/meschedl/Larval-Oyster-CASE-RNA/blob/master/scripts/Signed-Neg-Log-Pvalue.R) script if you want to know how to get those) for ALL the genes. And they need to be separated by a comma, or a csv file.

Cut the 1st (gene name) and the 10th (-log signed pvalue) columns in the new results file,

```
cut -f1,10 results_CA_F_log.txt > res_CA_log.txt
cut -f1,10 results_CASE_F_log.txt > res_CASE_log
cut -f1,10 results_SE_F_log.txt > res_SE_log
```

The files for some reason have "|LOC*" where the "*" stands for some number. The STR.GO.list below does not have this. I will try to get rid of the "LOC" string with the code below. 

```
awk '{sub(/\|LOC[0-9]+/, "\t", $1)}1' res_CA_log.txt > results_noLOC.txt
```

Then replace the tab with a comma.
```
sed -e 's/\t/,/g' results_noLOC.txt > res_CA_comma.csv
sed 's/,\s*/,/g' res_CA_comma.csv > res_CA_comma_nospace.csv

sed -e 's/\t/,/g' res_CASE_log > res_CASE_comma.csv
sed -e 's/\t/,/g' res_SE_log > res_SE_comma.csv
```


Then it needs a table with the first column being the gene name, and the second being the GO annotations associated with the gene. If there isn't an associated GO term it has to say unknown. We have a file with all the GO terms associated with the gene names, and even MSTRG transcripts ([J. Puritz](https://github.com/jpuritz) is a genius), but there are just blanks for no known GO term. 

The `STR.GO.list` file is located in `/home/Genomic_Resources/C_virginica/STR.GO.list`.

To make all the blank spaces unknowns:
```
mawk '{if (length($2) == 0) print $1 "\tunknown";
else print $0}' STR.GO.list > STR.GO.list.unknowns
```

Check for and remove a weird leading semicolon from in front of some of the GO terms
```
grep -o $'\t;GO' STR.GO.list.unknowns | wc -l
sed 's/\t;GO/\tGO/g' STR.GO.list.unknowns > STR.GO.list.unknowns.rmsemicolon
grep -o $'\t;GO' STR.GO.list.unknowns.rmsemicolon | wc -l
```

Then to get the GO term for all the ones in each results file, first cut the 1st column (gene name) from the results, then search that list in the GO list and make that output (the lines that match the searched for terms) to a new file for GO_MWU input.



```
cut -f1 results_noLOC.txt > res_CA_gene
grep -wFf res_CA_gene STR.GO.list.unknowns.rmsemicolon > res_CA_GO


cut -f1 res_SE_log > res_SE_gene
grep -wFf res_SE_gene STR.GO.list.unknowns.rmsemicolon > res_SE_GO
cut -f1 res_CASE_log > res_CASE_gene
grep -wFf res_CASE_gene STR.GO.list.unknowns.rmsemicolon > res_CASE_GO
```