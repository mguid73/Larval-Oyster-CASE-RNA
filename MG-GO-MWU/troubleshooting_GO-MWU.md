# Issues with GO MWU script

When I ran the `gomwuStats()` function on each of the datasets (CA,SE, and CASE) I received similar output.

*Note where it says at the bottom that 0 GO terms passed 10% FDR. The GO-MWU.R script says not to continue on if you get this message...*
But I feel like this isn't the case because Maggie had done this analysis previously and it did work so maybe there is something wrong with the input files?? 
Will dig into that next...

```
go.obo res_CA_GO res_CA_comma.csv BP largest=0.1 smallest=5 cutHeight=0.25

Run parameters:

largest GO category as fraction of all genes (largest)  : 0.1
         smallest GO category as # of genes (smallest)  : 5
                clustering threshold (clusterCutHeight) : 0.25

-----------------
retrieving GO hierarchy, reformatting data...

-------------
go_reformat:
Genes with GO annotations, but not listed in measure table: 0

Terms without defined level (old ontology?..): 2
-------------
-------------
go_nrify:
206 categories, 4401 genes; size range 5-2200.5
	3 too broad
	37 too small
	166 remaining

removing redundancy:

calculating GO term similarities based on shared genes...
72 non-redundant GO categories of good size
-------------

Secondary clustering:
calculating similarities....
Continuous measure of interest: will perform MWU test
**0 GO terms at 10% FDR**
```


For some reason the 'res_CA_GO' file (the annotation file) sometimes has semicolons after the tab for some of the GO terms. There are supposed to be semicolons separating mulitple GO terms per gene or MSTRG but I dont think that the inital GO term should have a leading semicolon. 

Example:
```
MSTRG.10004	;GO:0005576;GO:0005856;GO:0008092;GO:0008150;GO:0032991
MSTRG.10008	;GO:0003674;GO:0007165
MSTRG.10009	;GO:0003674;GO:0007165
MSTRG.10012	GO:0008233
...
```

I think I'll re-run the script that makes these files and see if that fixes it?

I got rid of the leading semicolon but that did not fix it...


NOW, I think the issue is that the MSTRG/gene names are not matching up between the .csv file and the GO annotation file. For instance, this annotation, "MSTRG.37002|LOC111106954" appears in column 1 of the .csv file but that annotation does not appear in the GO annotation (tab file). 