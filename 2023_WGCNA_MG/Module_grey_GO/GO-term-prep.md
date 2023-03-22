Find the genes in module grey that have "LOC" numbers associated with them using the following:

Use the csv of the genes in module grey as the input file and grep the .gtf file for LOC#

```
grep -f "module_grey_genes.csv" C_vir_merged.gtf | grep -oh "LOC[0-9]*" > greyloc

# get the unique occurences of those LOCs and store them in a csv
uniq greyloc >> unique_greyloc.csv
```

This provides a list of LOCs that are in module grey. We can then put these LOCs into gProfiler or GO-MWU.
