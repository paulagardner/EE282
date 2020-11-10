# Homework 3
wget chromosomes: wget, pasted the url from the search bar, and after /fasta: type /(specific filename we want) to pick files from the folder we're downloading from
```
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/dmel-all-chromosome-r6.36.fasta.gz

wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/md5sum.txt
```
md5 check: 
```
md5sum -c md5sum.txt

```
Now gunzip the .gz chromosome file we've selected to work with 
```
gunzip dmel-all-chromosome-r6.36.fasta.gz
faSize dmel-all-chromosome-r6.36.fasta

143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
*bases, sequnces, etc self-explanatory. N is an unknown nucleotide, and U and L are naming conventions for things you probably won't need anytime soon.*


now repeat wget for the gtf assembly files. In this case, we can use the asterisk to specify that we want all matches to this path, which we do, since the gtf folder only contains md5sum and the .gtf file. 
```
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/gtf/*
```
Note that here, use md5sum.txt.**1** since a md5sum.txt already exists 
```
md5sum -c md5sum.txt.1
```
formatted along Lauren's solutions: filtering the gtf file by the third column (where the features are), and making a new file from that
```
gawk -F '\t' '{print $3}' dmel-all-r6.36.gtf > col3_chrom.txt
```
Sort
```
sort col3_chrom.txt > col3_sorted.txt
```
get counts of the sorted column
```
uniq -c col3_sorted.txt > uniq_col3.txt
```
sort the counts in descending order
```
sort -k1,1nr uniq_col3.txt >uniq_col3_sorted.txt
```
Or, do it all as one line: 

```
gawk -F '\t' '{print $3}' dmel-all-r6.36.gtf | sort | uniq -c | sort -k1,1nr > piped.uniq.sort3.txt
```


Paavan's solution for number 2: more concise, but less exact - using word boundaries (?) this way. 

```
gawk -F '\t' '$3=="gene"{print $1}' dmel-all-r6.36.gtf | grep -P '^[4XY23][LR]?$' | sort | uniq -c
```
*Notes:  -F '\t': specifying what's happening with regards to the tabular nature of the data.*

*'$3 =="gene"{print$1}' filename ;  $ indicates columns, so this tells bash to filter column 1 (chromosome name) by the results of column 3 entries containing the word gene*

*grep -P '^[4XY23][LR]?$' takes that subset and filters for the first letter being an X,Y,2,3, etc. Second letter L or R*





# Notes 
JJ's solution for number 2: more specific, less concise- really spelling it out with the word boundaries front 
```
gawk -F '\t' ' $3 == 
"gene" {print $1 } ' <(zcat ee282/dmel-all-r6.36.gtf.gz) |sort | grep -P '^ (X|Y|(2L)|(2R)|(3L)|(3R)|4)$'| uniq -c | less
```

Starting out: 
All I did was navigate to myrepos/ee282/ , check what was going on with gitstatus, switch the branch I 
was working on to homework3 (which I had created earlier) and I was set!

```
[pdgardne@login-i17:~/myrepos/ee282] $ wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
```


Interestingly, although I expected gitstatus to tell me this download was on the homework3 branch; 
instead it went to the directory itself. I'm still fuzzy on how it all works:

```(base) [pdgardne@login-i17:~/myrepos/ee282] $ git status
#On branch homework3
nothing to commit, working directory clean
(base) [pdgardne@login-i17:~/myrepos/ee282] $ ls
code  data  dmel-all-chromosome-r6.36.fasta.gz  output  README.md  tmp
```

Next it's checksum time. Flybase, according to lecture 6, uses md5 as checksums. 
https://www.howtoforge.com/linux-md5sum-command/


```
[pdgardne@login-i17:~/myrepos/ee282] $ md5sum dmel-all-chromosome-r6.36.fasta.gz
fbd2855a20c3610050ff2496dd975821  dmel-all-chromosome-r6.36.fasta.gz
```

```
[pdgardne@login-i17:~/myrepos/ee282] $ md5sum -c dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.36_FB2020_05/fasta/md5sum.txt
```
https://www.oreilly.com/library/view/linux-pocket-guide/9780596806347/re69.html. 

```
md5sum -c md5sum.txt shows all the checksums that match
md5sum *.gz > md5check.txt
diff md5{sum,check}.txt
```

#you get a different output than JJ did in his tutorial, but note that he used wget to get ALL the fastas in 
#that file location: wget -R "*aligned*" \ ftp://ftp.flybase.net/releases/FB2020_05/dmel_r6.36/fasta/*


#NOW: Calculate Summaries of the Genome
#Total number of nucleotides
#Total number of Ns
#Total number of sequences
#was given a hint for the assignment in general to use fasize. "Also, investigate out the gtf file format specification"
#https://uswest.ensembl.org/info/website/upload/gff.html

don't forget to unzip!
gunzip dmel-all-chromosome-r6.36.fasta.gz

Ok update: for when you get to this part, the 
internet will tell you that 'grep ">" dmel-all-chromosome-r6.36.fasta' will work- it does NOT

[pdgardne@hpc3-19-00:~/myrepos/ee282] $ grep -v ">" dmel-all-chromosome-r6.36.fasta | wc
1797519 1797496 145523521

#Question: Useful Bash Commands To Handle Fasta Files
https://www.biostars.org/p/17680/ , which recommends awk.
awk help

