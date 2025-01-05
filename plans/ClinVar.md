## ClinVar

Some of the actionable findings should be restricted to phenotype-specific pathogenic variants. This is a little tricky, because
while a clinvar entry (allele ID) may have multiple associated phenotypes, I lost this in the aggregation process.

We can recover some of this information by pulling the variant.txt file from the ClinVar FTP site, then filter for a combination of disease name and gene ID/Symbol. This will deliver a collection of AlleleIDs, or even chrom/position/ref/alt, either of which can be used to filter the already-applied ClinVar annotations.

This seems easier than renovating the ClinvArbitration process to do some new bespoke things.

This is an inherently more complicated decision tree than the A (dominant) and B (recessive) variant types, so saving that for after an initial MVP.
