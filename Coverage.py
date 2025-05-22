#Starting with the program of slide 19 and working with BAM file
#10_Normal_Chr21.bam,  find the locus with:
#   Single diploid-organism heterozygosity (approx.), and highest coverage
#   Output the locus, its alleles, and the number of high-quality reads for each
#   allele.
#Then, add the code to filter out bad/poor alignments (slide 20) to your program.
#   Does the highest coverage heterozygous locus and its read counts change? If
#   so, how?

import pysam
bf = pysam.Samfile('10_Normal_Chr21.bam')

max_coverage = 0                #stores highest coverage found
max_heterozygosity_pos = None   #stores position of locus with highest heterozygosity
max_counts = {}                 # stores counts of alleles in max heterozygosity position

for pileup in bf.pileup('21'):
    counts = {}
    # ...examine every aligned read
    for pileupread in pileup.pileups:
        # ...check the read and alignment
        if pileupread.indel:
            continue
        if pileupread.is_del:
            continue
        al = pileupread.alignment
        if al.is_unmapped:
            continue
        if al.is_secondary:
            continue
        if int(al.opt('NM')) > 1:
            continue
        if int(al.opt('NH')) > 1:
            continue
        # ...and get the read-base
        if not pileupread.query_position:
            continue
        readbase = al.seq[pileupread.query_position]
        if not pileupread.query_position:
            continue
        readbase = pileupread.alignment.seq[pileupread.query_position]
        # Count the number of each base
        if readbase not in counts:
            counts[readbase] = 0
        counts[readbase] += 1
    # If there is no variation, move on
    if len(counts) != 2:                            #skips if there are not exactly 2 different bases to indicate heterozygosity
        continue

    #threshold difference
    allele_count = list(counts.values())            #list of counts of each allele
    diff = abs(allele_count[0] - allele_count[1])   #absolute different between counts of two alleles
    avg = (allele_count[0] + allele_count[1])/2     #average of the two counts
    percent_avg = (diff/avg)*100                    #% difference between the two allele counts
    
    if percent_avg < 20:                            #checks if the percent average is less than 20%
        #coverage
        coverage = pileup.n                         #gets the coverage (# of reads) at the pileup position
        if coverage > max_coverage:                 #checks if its greater than max coverage
            max_coverage = coverage                 #updates these values to that of the highest coverage locus found
            max_heterozygosity_pos = pileup.pos
            max_counts = counts
            #print(pileup.pos, pileup.n, end = " ")
            #for base in sorted(counts):
            #    print(base,counts[base], end = " ")
            #print()

#output position with max heterozygosity and its coverage
if max_heterozygosity_pos is not None:
    print("Locus Position:", max_heterozygosity_pos)
    print("Max Coverage:", max_coverage)
    for alleles, counts in sorted(max_counts.items()):
        print("Allele:", alleles, "Reads:", counts)
