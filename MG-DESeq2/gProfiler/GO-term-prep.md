Find the DEGs in each dataset (CA, SE, CASE) that have "LOC" numbers associated with them using the following:

Use the csv of the significant, ordered genes (from most to least differentially expressed by absolute value of the log2FoldChange) as the input file and grep the .gtf file for LOC#

```
# first extract the first column from "ordered_lfc_CA.txt" to use as the input file for the grep 
cut -f1 ordered_lfc_CA.txt > CA_DEGs.txt

# grep that file for the LOCs and save as a file
grep -f "CA_DEGs.txt"  ../gProfiler/C_vir_merged.gtf | grep -oh "LOC[0-9]*" > CAloc

# get the unique occurences of those LOCs and store them in a csv
uniq CAloc >> unique_CAloc.csv
```

Use the same code to prep the LOCs for all treatments.

I ended up with 3 files (1 per treatment):

```
unique_CAloc.csv
unique_SEloc.csv
unique_CASEloc.csv
```

This provides a list of LOCs that are associated with each treatment. We can then put these LOCs into gProfiler or GO-MWU.

[gProfiler for CA DEGs](https://biit.cs.ut.ee/gplink/l/HJa14VFPTN)
- **detection of a chemical stimulus (GO:0009593)**
- **sensory perception of a chemical stimulus (GO:0007606)**
- **G-protein coupled receptor activity**

[gProfiler for SE DEGs](https://biit.cs.ut.ee/gplink/l/ji8YJMGZTU)
- **metal cluster binding (GO:0051540)**: "Binding to a cluster of atoms including both metal ions and nonmetal atoms, usually sulfur and oxygen. Examples include iron-sulfur clusters and nickel-iron-sulfur clusters."
- **metal ion binding (GO:0046872)**
- **organonitrogen compound metabolic processes (GO:1901564)**: "The chemical reactions and pathways involving organonitrogen compound."
- **ethanol metabolic process (GO:0006067)**
- **ethanol oxidation (GO:0006069)**
- **sulfate assimilation (GO:0000103):** sulfate assimilation, phosphoadenylyl sulfate reduction by an oxidoreductase, acting on sulfur group of donors, NAD or NADP as acceptor | sulphate assimilation | sulphate assimilation, phosphoadenylyl sulphate reduction by an oxidoreductase, acting on sulphur group of donors, NAD or NADP as acceptor
- **water soluble vitamin biosynthetic processes**


[gProfiler for CASE DEGs](https://biit.cs.ut.ee/gplink/l/jFl3OYoTTs)
- **oxidoreductase activity (GO:0016491):** Catalysis of an oxidation-reduction (redox) reaction, a reversible chemical reaction in which the oxidation state of an atom or atoms within a molecule is altered. One substrate acts as a hydrogen or electron donor and becomes oxidized, while the other acts as hydrogen or electron acceptor and becomes reduced.
- **phenylalanine 4-monooxygenase activity (GO:0004505):** Catalysis of the reaction: L-phenylalanine + tetrahydrobiopterin + O2 = L-tyrosine + 4-alpha-hydroxytetrahydrobiopterin. Phenylalanine is an essential amino acid for humans (meaning that we must get it from teh foods we eat). Lots of suger-substituted foods contain phenylalanine. 
- **metal cluster binding GO:0051540**
- **iron-sulfur cluster binding GO:0051536**
- **L-phenylalanine catabolic process GO:0006559 GO:0006558**
- **erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process GO:1902221**: these are precursor molecules to phenylalanine
- **cholesterol transport/sterol transport**
- **PcG protein complex (GO:0031519):** A chromatin-associated multiprotein complex containing Polycomb Group proteins. In Drosophila, Polycomb group proteins are involved in the long-term maintenance of gene repression, and PcG protein complexes associate with Polycomb group response elements (PREs) in target genes to regulate higher-order chromatin structure. PMID:9372908


## But which terms are shared between treatments and which are unique to those treatments?

To compare GO terms between the groups, let's look at the unique_'treatment'loc.csv files and see which GO terms appear in multiple datasets or are unique to their given dataset.

First let's find the LOCs that appear in all 3 files.
```
cd ../outputs
grep -xF -f unique_CAloc.csv unique_SEloc.csv | grep -xF -f - unique_CASEloc.csv > sharedbyall

uniq sharedbyall > sharedbyall.csv
```

## GO terms shared between CA and SE but not in CASE
(the grep -v flag indicates 'not shared')
```
grep -xF -f unique_CAloc.csv unique_SEloc.csv | uniq > sharedby_CA_SE
grep -xFv -f unique_CASEloc.csv sharedby_CA_SE > sharedby_CA_SE_uniq.csv
```
There are seven LOCs shared by CA and SE that aren't also in CASE.


## GO terms shared by CA and CASE but not in SE
```
grep -xF -f unique_CAloc.csv unique_CASEloc.csv | uniq > sharedby_CA_CASE
grep -xFv -f unique_SEloc.csv sharedby_CA_CASE > sharedby_CA_CASE_uniq.csv
```
There are 182 LOCs shared by CA and CASE that arent also in SE

gProfiler output: `gProfiler_shared_CAandCASE.csv` 
Nothing really significant or super interesting in here at first glance


## GO terms shared by SE and CASE but not in CA
```
grep -xF -f unique_SEloc.csv unique_CASEloc.csv | uniq > sharedby_SE_CASE
grep -xFv -f unique_CAloc.csv sharedby_SE_CASE > sharedby_SE_CASE_uniq.csv
```
There are 43 LOCs shared by SE and CASE that arent also in CA
gProfiler output: `gProfiler_shared_SEandCASE.csv`


## GO terms in CASE but not in CA and SE
```
grep -xFv -f sharedby_CA_SE unique_CASEloc.csv > unique_to_CASE.csv
```
355 LOCs in CASE but not in CA and SE

**Summary: There are GO terms associated with electron transfer, phenylalanine metabolism, steriod/hormone metabolism/activity**

terms of interest from `gProfiler_uniquetoCASE_output.csv`

- "GO:MF","oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced pteridine as one donor, and incorporation of one atom of oxygen","GO:0016714"
- "GO:MF","phenylalanine 4-monooxygenase activity","GO:0004505"
- "GO:MF","nuclear steroid receptor activity","GO:0003707"
- "GO:MF","nuclear steroid receptor activity","GO:0003707"
- "GO:MF","GTPase activator activity","GO:0005096"
- "GO:MF","oxidoreductase activity","GO:0016491"
- "GO:BP","small molecule metabolic process","GO:0044281"
- "GO:BP","L-phenylalanine catabolic process","GO:0006559"
- "GO:BP","erythrose 4-phosphate/phosphoenolpyruvate family amino acid metabolic process","GO:1902221"
- "GO:BP","L-phenylalanine metabolic process","GO:0006558"
- "GO:BP","cholesterol transport","GO:0030301"
- "GO:BP","intracellular cholesterol transport"
- "GO:BP","enteroendocrine cell differentiation","GO:0035883"
- "GO:BP","endocrine system development","GO:0035270"
- "GO:BP","peptide hormone secretion","GO:0030072"