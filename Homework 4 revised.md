# Homework 4
get flybase genome: 
```
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/dmel-all-chromosome-r6.36.fasta.gz
```



*beware! you spent like an hour here with some weird typos that were causing errors; mostly "command not found"- mostly just honest typos and some possibly whitespace issues?*



**partition the genome into sequences greater than and less than/equal to 100000 bp:**

**Less than:**
```
bioawk -c fastx \
' length($seq) <= 100000 { print ">" $name "\n" $seq } ' \
dmel-all-chromosome-r6.36.fasta.gz \
| gzip -c \
> dmelshortsequences.fasta.gz

faSize dmelshortsequences.fasta.gz
```

**Greater than:**
```
bioawk -c fastx \
' length($seq) > 100000 { print ">" $name "\n" $seq } ' \
dmel-all-chromosome-r6.36.fasta.gz \
| gzip -c \
> dmellongsequences.fasta.gz

faSize dmelshortsequences.fasta.gz
```
*reference your unedited version for an FaFilter way to do this*



## **Calculate 3 metrics for all sequences-- FaSize will return them** 

*first demos: length/gc on their own
partitioning for sequence lengths and gc contents* 

**Making new partitions:** *(could pipe from previous by deleting the condition filter (lenth($seq) <= 100000), and piping from the already partitioned file, but anyway)*:

* Short sequences
```
bioawk -c fastx \
' length($seq) <= 100000 { print ">" $name "\n" length($seq) "\t" gc($seq) } ' \
dmel-all-chromosome-r6.36.fasta.gz \
| gzip -c \
> dmelshortsequences.info.fasta.gz \
|faSize dmelshortsequences.info.fasta.gz
```
* Large sequences
```
bioawk -c fastx \
' length($seq) > 100000 { print ">" $name "\n" length($seq) "\t" gc($seq) } ' \
dmel-all-chromosome-r6.36.fasta.gz \
| gzip -c \
> dmellongsequences.info.fasta.gz \
|faSize dmellongsequences.info.fasta.gz
```


# **PLOTS**

*following the demo of doing it as a pipe from the our partitioned files*
*reference unedited version of this file here to see lots of extraneous code (from demos) for CDF plots that tell you about your data/what's going on with the files, but we instead want histograms from R*


* *need to save the sorted data as CSV so R can read them. You did the same thing where you made new partitions again, but remember; you can pipe from the original partition files easily too*
  
**Less than:**
 ```
bioawk -c fastx \
' length($seq) <= 100000 { print length($seq) "," gc($seq) } ' \
dmel-all-chromosome-r6.36.fasta.gz \
> dmelshortsequences.csv
```


**Greater than:**

 ```
bioawk -c fastx \
' length($seq) > 100000 { print length($seq) "," gc($seq) } ' \
dmel-all-chromosome-r6.36.fasta.gz \
> dmellongsequences.fasta.csv
```
Use webdrive to get the .csv files to R (just copied them into a folder on your own hard drive)

**R CODE:**
* *see unedited version of this file for R code with notes, about how to make nicely formatted histograms- however, the very straightforward way to do this, since we're interested in quick data vis. Dplyr, tibbles vs data.frames, etc*
```
hist(log(dmelshortsequences[[1]]))
#got log() idea from Karina. Thanks!! This is a graph of sequence length distribution for the less than 100000 bp file

hist(dmelshortsequences[[2]])
#short sequences gc% distribution

hist(dmellongsequences[[1]]) 
#no log needed here, since there are far fewer long sequences. Length distribution

hist(dmellongsequences[[2]])
#gc % distribution for long sequences
```


## **Cumulative plotCDF**


**Less than:**
* *start sorting on column 1, end sorting on column 1, reverse sort with r, turn on numerical sorting by doing n*
```
bioawk -c fastx \
' { print length($seq) "\t" gc($seq) } ' \
dmelshortsequences.fasta.gz \
| sort -k1,1rn \
> dmelshortsequences.lengthgc.lengthsort.txt
```
* Plot of cumulative sequence lengths (of short sequences) when the file was sorted largest to smallest: 
```
plotCDF <(cut -f 1 dmelshortsequences.lengthgc.lengthsort.txt) dmelshortsequences.lengthgc.lengthsort.png

display dmelshortsequences.lengthgc.lengthsort.png
```


**Greater than:**
```
bioawk -c fastx \
' { print length($seq) "\t" gc($seq) } ' \
dmellongsequences.fasta.gz \
| sort -k1,1rn \
> dmellongsequences.lengthgc.lengthsort.txt
```

* Plot of cumulative sequence lengths (of short sequences) when the file was sorted largest to smallest: 

```
plotCDF <(cut -f 1 dmellongsequences.lengthgc.lengthsort.txt) dmellongsequences.lengthgc.lengthsort.png\
display dmellongsequences.lengthgc.lengthsort.png
```


# **Genome Assembly**
*in unedited file, lots of instructions that are for following along in a demo where JJ renamed files and specified paths as variables(?), but here's the most straightforward way to do it*
```
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz
```
* *interestingly, it comes up as this file having been created in 2018- from the directory I got it from, presumably* 
```
conda activate ee282
conda install -y miniasm minimap
```
**get 32 cores for your work- remember to log out when you're done!!** (with exit)
```
srun -c 32 -A ecoevo282 --pty --x11 bash -i
```

* *this helps you work w/ nanopore reads- the -Sw5, -L100, -m0 command is a preset that works well for doing oxford nanopore*
```
minimap -t 32 -Sw5 -L100 -m0 iso1_onp_a2_1kb.fastq.gz{,} \
| gzip -1 \
> processed.gz
```
*the {,} command is important in bash- we can use brace expansion here to map the read against itself- target and query are the same. simplest way to tell bash to do that*

* *miniasm figures out what to do with overlaps*
```
miniasm -f iso1_onp_a2_1kb.fastq.gz processed.gz \
> processed.gfa
```


## **Assembly assesment**
* N50 number: how long is the sequence at a certain point in the graph-- proxy for how contiguous the genome is

```
n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
} 

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' processed.gfa \
| tee >(n50 /dev/stdin > n50.txt) \
| fold -w 60 \
> processed.fa
```

* *awk command: for every line that has an S at the beginning, print (in fasta format), print $2, make a new line, then print $3 (which are name and then sequence)*
* *tee lets you turn reads.gfa into the reports file*
* *fold takes very  long sequences (not having lines millions of sequences long), and will keep the lines at 60 characters long*
* N50 demos are in nov 30 and dec 4 class videos


### **Compare your N50 against that of the drosophila community's:**
checking the result:
```
faSize -detailed $procesed/processed.fa | less
```

so the N50 score for the assembly is a lot larger than the one I made- which makes sense, 
since this is kind of the 'quick and dirty' way to do an assembly- it follows that the minimum size is smaller for my assembly than the standard one. 




### **CDF plot**
* FAsplitbyN takes your scaffolds and turns it into contigs
  * lets you orient all the contig constituents of a scaffold, but faSplitbyN splits them up so you can directly compare your contig to the reference one

* after Karina's code:
```
bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' miniasm_assembly_N50processed.fa \
| sort -k1,1rn \
> miniasm_assembly_length.txt

plotCDF <(cut -f 1 miniasm_assembly_length.txt) miniasm_assembly_contigplot.png

display miniasm_assembly_contigplot.png
plotCDF <(cut -f 1 assembly_length.txt) assembly_contigplot.png \
display assembly_contigplot.png


bioawk -c fastx \
  ' { print length($seq) "\t" gc($seq) } ' miniasm_assembly_N50processed.fa \
| sort -k1,1rn \
> miniasm_assembly_length.txt

plotCDF <(cut -f 1 miniasm_assembly_length.txt) miniasm_assembly_contigplot.png

display miniasm_assembly_contigplot.png
plotCDF <(cut -f 1 assembly_length.txt) assembly_contigplot.png \
display assembly_contigplot.png


```

# BUSCO analysis- following tutorial from university of Geneva
*busco can be a good way to assess completeness of your assembly*
*jj searches bioconda busco recipe*
conda install busco

*busco assumes genes it's looking for are always present, and only present in one copy-- no gene duplication- so obviously must interpret results with this in mind*
*orthologs- genes shared through common descent. Contrast w/ paralog- so you may want to compare against lineages*

busco -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]