# Coverage-for-Single-Diploid-Organism-Heterozygosity

The highest coverage heterozygous locus and its read counts changed after the bad/poor alignments had been filtered out.  The locus position changed, the max coverage decreased, and the reads both increased.  The change occurred as:

Without Slide 20 Code:
Locus Position: 33040370
Max Coverage: 1050
Allele: C Reads: 2
Allele: T Reads: 2

With Slide 20 Codes: 
Locus Position: 46270195 
Max Coverage: 747
Allele: A Reads: 316
Allele: G Reads: 350

The results change due to the filter skipping over indels, reads marked as deleted, unmapped reads, secondary alignments, mismatches, reads mapping to multiple locations, and reads whose query position is invalid.  The code filters out reads that are likely to introduce noise or inaccuracies into the analysis. This enhances the reliability of the results when determining heterozygosity and counts at specific loci.
