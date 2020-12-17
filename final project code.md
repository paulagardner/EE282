# METHODS
In order to compare the gc content per chromosome per species, and thus have some idea of how gc content compares to genome size in different clades, and how the clades stack up against each other, it was necessary to get assemblies that would stand in as representative data points for different clades. I took the assemblies of five different species: an invertebrate, *Drosophila melanogaster*, a tetrapod outgroup, the fish *Danio rerio*, a bird, to stand in for sauropsids, *Taeniopygia guttata*, and two mammals: a marsupial, *Monodelphis domestica*, and the eutherian mammal *Mus musculus*. I wrote various lines of code and pipelines to look at and then partition the data in the bash terminal. First, I looked at which fasta lines contained the most complete contig for any given chromosome, and used bioawk to filter the chromosomes based on that, giving me .csv files that contained three columns: some form of chromosome ID, the sequence length for the chromosome contig, and the gc proportion (expressed as a decimal). In the case of the zebra finch, *T. gutatta*, it was necessary because of less complete chromosome contigs to filter by name, as in this case the longest chromosomes were not always the most representative ones. 

The .csv files were then imported into R studio, where they were collected into one large data frame with an associated ID, which allowed for a color-coded scatterplot of length (along the x axis) and gc percentage (along the y axis) to be created. Following the work of Li and Xiu-Qing, A test of Spearman's correlation coefficient was applied to compare the length and gc content among chromosomes. 

# RESULTS
## Results for spearman's correlation coefficient: 
```
data:  fly$length and fly$`gc%`
S = 14, p-value = 0.2417
alternative hypothesis: true rho is not equal to 0
sample estimates:
rho 
0.6 


	
data:  fish$length and fish$`gc%`
S = 2008, p-value = 0.5529
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.1269565 


data:  zebrafinch$length and zebrafinch$`gc%`
S = 13094, p-value = 1.375e-07
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.8338936 


data:  opossum$length and opossum$`gc%`
S = 94, p-value = 0.793
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.1190476 


data:  mouse$length and mouse$`gc%`
S = 1576, p-value = 0.4333
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.1849624 
```


[<img src="https://github.com/paulagardner/EE282_back/blob/final_project/scatterplot-gc%25vslength.png">](http://github.com)
See above or in the github repository for the scatterplot of length plotted against GC content. Each point represents an individual chromosome.  There seems to be little correlation between chromosome length and GC content in these animals. 
Most instances of computing spearman's correlation coefficient for the species sampled  resulted in p values > .05


# DISCUSSION 
Broadly speaking, I was interested in how following the analysis of the Li and Xiu-Qing paper might translate to looking at an actual dataset. This is my first instance of actually trying to replicate a set of methods I took from a paper, and the results are necessarily limited and are otherwise generally smaller in scope. 

The only significant finding for spearman's correlation coefficient within my samples was a strongly negative relationship between gc content as chromosome size went up in the Zebra finch (* T. guttata*), which is in line with the Xiu et. al paper. Interestingly, although the other species had very large p-values, this negative correlation held true for the 'higher order' vertebrates; and was not true for the fish and fly assemblies. Also, it's worth noting that the zebra finch had the largest number of chromosomes. It would be interesting to see if comparing two different sets of assemblies with the same range of phyla but with smaller or larger chromosome numbers would have much of an effect. 

The authors speak a lot about the data in terms of averages, so my one-sample-per-phyla approach means direct comparison is not possible. However, I am skeptical of many of the claims, since an 'average' is so highly dependent on how the phyla would be sampled. The authors state that "C+G content...was highly similar among species of mammalian animals...and also clearly different between non-mammalian animals and mammalian animals". Although they were speaking about averages in the context of total genomic length, this is something that doesn't seem to hold true for my data, as the marsupial's GC content per chromosome (which is still a measure of GC content over length) was markedly different from the mouse's. Their data does not include any non-eutherian mammals, and I take umbrage with their discussing average genome length for mammals vs. non-mammals when nearly half of the eutherians sampled were primates. 

The scatterplot I produced should be taken with a grain of salt; since in some cases a lot of non-contiguous genomic data was discarded. It's interesting to note anyway that, at least in my sample, there was little correlation between how derived the organism was on the tree of life vs. its chromosome number, or otherwise its chromosome length. The marsupial had far and away the largest chromosomes and a middling GC content within them, versus its next closest relative, the mouse, having overall the lowest GC content per chromosome. If GC content scaled proportionally with chromosome length, we would expect a line bisecting the graph from the origin to its upper right corner- which is clearly not the case. Overall, there is little correlation between the 'derived-ness' of a phylum to which a sampled organism belonged and the characteristics of its chromosomes in terms of GC proportion (used as a proxy for gene location/proportion) and length. 


# CODE BELOW: 

## *Danio rerio*, zebrafish

* get the assembly: 
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/002/035/GCA_000002035.4_GRCz11/GCA_000002035.4_GRCz11_genomic.fna.gz 
```
* Below is code from JJ- to look at which chromosomes exist on this file, and how long each are so as to be able to only pull out the "real" chromosomes. In the case of this danio rerio file, that required looking at the NCBI website to see which codes corresponded to chromosome contigs, as they weren't labed chr_, for example. 
```
faSize -detailed GCA_000002035.4_GRCz11_genomic.fna.gz \
| sort -k2,2rn \
| less -N
```
* partition sequences greater than 1000000 base pairs and print out three columns: name, sequence length, and gc%. Sort in descending order.
```
bioawk -c fastx \
' length($seq) > 1000000 { print $name "\t" length($seq) "\t" gc($seq) } ' \
GCA_000002035.4_GRCz11_genomic.fna.gz \
|sort -k1,1nr \
| gzip -c \
> danrerchromosomes.info.fasta.gz 
```
* same partitioning, but instead of saving as a .gz, save as a comma separated, csv. (JJ has since reminded me that R can read tabulated data, so this part was unnecessary, but nonetheless I did it for all species)
```
bioawk -c fastx \
' length($seq) > 1000000 { print ($name) "," length($seq) "," gc($seq) } ' \
GCA_000002035.4_GRCz11_genomic.fna.gz \
|sort -k1,1nr \
> danrerchromosomes.info.csv 

```

## *Drosophila melanogaster* (fruit fly)
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.masked.gz

faSize -detailed  dm6.fa.masked.gz \
| sort -k2,2rn \
| less -N


bioawk -c fastx \
' length($seq) > 1000000 { print $name "\t" length($seq) "\t" gc($seq) } ' \
dm6.fa.masked.gz \
| gzip -c \
> dmelchromosomes.info.fasta.gz 

bioawk -c fastx \
' length($seq) > 1000000 { print ($name) "," length($seq) "," gc($seq) } ' \
dm6.fa.masked.gz \
|sort -k1,1nr \
> dmelchromosomes.info.csv 

```





## *Taeniopygia guttata*, zebra finch

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/taeGut2/bigZips/taeGut2.fa.gz
 
faSize -detailed  taeGut2.fa.gz \
| sort -k2,2rn \
| less -N
```
#reviewing this, and there are some files that are not 'chromosome' contigs that are bigger than some
#of the other chromosomes you want to match to. bioawk then grep? 

#need to figure out regular expressions

```
VERSIONS: 
(ee282) [pdgardne@hpc3-17-11:~/myrepos/ee282/finalproject] $ perl -lane ' if ($_ =~ /^([^_\t]+)_?([^\t]*)\t(.*)$/) { print $1, "\t", $2, "\t", $3 }  '   /dfs6/pub/jje/ee282/tmp.txt |less

(ee282) [pdgardne@hpc3-17-11:~/myrepos/ee282/finalproject] $ perl -lane ' if ($_ =~ /^([^_\t]+)_?([^\t]*)\t(\d+)$/) { print $1, "\t", $2, "\t", $3 }  '   /dfs6/pub/jje/ee282/tmp.txt |less

(ee282) [pdgardne@hpc3-17-11:~/myrepos/ee282/finalproject] $ perl -lane ' if ($_ =~ /^([^_\t]+)_?([^\t]*)\t([0-9]+)$/) { print $1, "\t", $2, "\t", $3 }  '   /dfs6/pub/jje/ee282/tmp.txt |less

bioawk -c fastx ' { print $name, length($seq), gc($seq) } ' \
  /dfs6/pub/jje/ee282/taeGut2.fa.gz \
| sort -k2,2nr \
| perl -lane \
  ' if ($_ =~ /^([^_\t]+)_?([^\t]*)\t(.+)$/) { print $1, "\t", $2, "\t", $3 }  ' \
| less


#JJ's version: 
bioawk -c fastx ' { print $name, length($seq), gc($seq) } ' \
  taeGut2.fa.gz \
| sort -k2,2nr \
| perl -lane \
  ' if ($_ =~ /^([^_\t]+)_?([^\t]*)\t(.+)$/) { print $1, ",", $2, ",", $3 }  ' \
> taegutchromosomes.info.csv

bioawk -c fastx ' { print $name, length($seq), gc($seq) } ' \
  taeGut2.fa.gz \
| sort -k2,2nr \
| perl -lane \
  ' if ($_ =~ /^([^_\t]+)_?\t([^\t]*)\t(.+)$/) { print $1, ",", $2, ",", $3 }  ' \
> taegutchromosomes.info.csv
```


# MAMMAL: 

## *Mus musculus*, house mouse
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.masked.gz

faSize -detailed  mm39.fa.masked.gz \
| sort -k2,2rn \
| less -N


bioawk -c fastx \
' length($seq) > 1000000 { print $name "\t" length($seq) "\t" gc($seq) } ' \
mm39.fa.masked.gz \
|sort -k1,1nr \
| gzip -c \
> muschromosomes.info.fasta.gz 

bioawk -c fastx \
' length($seq) > 1000000 { print ($name) "," length($seq) "," gc($seq) } ' \
mm39.fa.masked.gz \
|sort -k1,1nr \
> muschromosomes.info.csv 
```


## *Monodelphis domestica*, Gray short-tailed opossum
```
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/monDom5/chromosomes/*

rm chrUn.fa.gz
rm chrM.fa.gz
```
* potentially shouldn't have removed chromosome M, but I did for my analysis. 
  
```
cat *fasta > opossum.fa.gz


faSize -detailed  opossum.fa.gz \
| sort -k2,2rn \
| less -N

cp opossum.fa.gz ~/myrepos/ee282/finalproject

bioawk -c fastx \
' { print ($name) "," length($seq) "," gc($seq) } ' \
opossum.fa.gz \
|sort -k1,1nr \
> possumchromosomes.info.csv 
```


R CODE: 

```
library(dplyr)

fly <- dmelchromosomes_info
colnames(fly) <- c("chr", "length", "gc%")

fish <- danrerchromosomes_info
colnames(fish) <- c("chr", "length", "gc%")

zebrafinch <- taegutchromosomes_info
colnames(zebrafinch) <- c("chr", "length", "gc%")

mouse <- muschromosomes_info
colnames(mouse) <- c("chr", "length", "gc%")

opossum <- possumchromosomes_info
colnames(opossum) <- c("chr", "length", "gc%")

library(ggplot2)

ggplot(aes(x='length', y="gc%"), data=fly) + 
  geom_point(aes(x='length', y= 'gc%'), data=fish)

#trying to add an ID column so I can make one big data frame: 

fly['ID'] = "Fly"
fish['ID'] = 'Fish'
zebrafinch['ID'] = 'Finch'
opossum['ID'] = 'Opossum'
mouse['ID'] = 'Mouse'

dataframe <- rbind (fish, fly, zebrafinch, opossum, mouse)

#this wasn't working, and the issue was the %. Using backticks allows you 
#to make names explicitly
plot <- ggplot(dataframe, aes(length, `gc%`, color=ID)) + geom_point()
plot
#plot + scale_x_continuous(trans='log2')


#after this link: http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
#spearman's correlation coefficient: 
install.packages("ggpubr")
library("ggpubr")

flycor <- cor.test(fly$length, fly$`gc%`, method=c("spearman"))
flycor

fishcor <- cor.test(fish$length, fish$`gc%`, method=c("spearman"))
fishcor

zebrafinchcor <- cor.test(zebrafinch$length, zebrafinch$`gc%`, method=c("spearman"))
zebrafinchcor

opossumcor <- cor.test(opossum$length, opossum$`gc%`, method=c("spearman"))
opossumcor

mousecor <- cor.test(mouse$length, mouse$`gc%`, method =c("spearman"))
mousecor
```

