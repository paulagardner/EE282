# Project Analysis proposal

As mentioned in my topic proposal, GC content in general and CpG islands in particular are a common metric used across clades as a metric for many different types of analysis, from having biochemical importance in the context of DNA methylation [1] to being indicators for gene locations in eukaryotes [2]. Exploring the possibilities of working with this type of data may prove useful due to this broad range of applications in the bioinformatics context.

The Li and Xiu-Qing paper [3] found a negative correlation between GC content and chromosome size in non-mammalian animals. I plan to reference sequences from the UCSC genome browser gateway, and look at 3 genome assemblies: the fruit fly (Drosophila Melanogaster) dm6 assembly; the Zebra finch (Taeniopygia guttata) taegut2 assembly; and the zebrafish (Danio rerio) danRer11 assembly. The referenced paper compared G+C content to average chromosome size, and claims a highly significant negative correlation in this clade between chromosome size and C+G content.

I plan to write a pipeline that will import the assemblies in the terminal, read through the sequences, and return the G+C content of each chromosome, and integrate that with other metrics- for instance, chromosome size. I would plan to try to make a pipeline myself for reading the sequence in each chromosome and calculating the proportion those G and C bases occupy of the total sequence, and take that output and compare it to chromosome length, and test that as the authors did with spearman's correlation coefficient.

The plan is to use a few conda packages to do the work, including fasize to calculate chromosome length and awk or bioawk to actually find the G+C counts. Time permitting, I'd add more reference assemblies and look into whether the negative correlation seems to grow stronger or weaker over the radiation from invertebrates to non-mammalian vertebrates.

If this all proves to be feasible to accomplish, then I might extend my analysis to include searching for CpG islands in some of the chromosomes of these assemblies by looking for presence of many C nucleotides immediately followed by a G nucleotides in selected chromosomes of the different species, and analyzing this content with the methods above.

## References

[1] Jang, Hyun Sik, Woo Jung Shin, Jeong Eon Lee, and Jeong Tae Do. 2017. “CpG and Non-CpG Methylation in Epigenetic Gene Regulation and Brain Function.” Genes 8 (6). <https://doi.org/10.3390/genes8060148>.

[2] Takai, Daiya, and Peter A. Jones. 2002. “Comprehensive Analysis of CpG Islands in Human Chromosomes 21 and 22.” Proceedings of the National Academy of Sciences of the United States of America 99 (6): 3740–45. <https://doi.org/10.1073/pnas.052410099>.

[3] Li, Xiu-Qing, and Donglei Du. 2014. “Variation, Evolution, and Correlation Analysis of C+G Content and Genome or Chromosome Size in Different Kingdoms and Phyla.” PloS One 9 (2): e88339. <https://doi.org/10.1371/journal.pone.0088339>.