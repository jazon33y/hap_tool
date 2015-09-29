hap_tool
========

This command-line tool (script) uses phylogenetic information to infer false negative rates in mitochondrial genome data from high-throughput sequencing experiments using estimates of haplogroup assignments using the haplogrep algorithm (doi: 10.1002/humu.21382).


```
usage: hap_tool.py [-h] [-ana <ana>] [-ft <ft>] -file [<file> [<file> ...]]

This script uses phylotree information to: 

	- return chrM in FASTA format
	- estimate haplogroup
	- estimate false negative call rates
	- return all loci which should have a SNP

optional arguments:
  -h, --help            show this help message and exit
  -ana <ana>            analysis type; fasta, haplogroup, FN_rate, FN_locus
  -ft <ft>              file type; vcf-hg19, vcf-grch37, var
  -file [<file> [<file> ...]]
                        input file, .vcf or .var
```

For example, you can use one of the tester files to estimate its haplogroup:

```
$ python hap_tool.py -ana haplogroup -ft vcf-hg19 -file FreeBayes.vcf
FreeBayes.vcf: K1a4a1a+195; 0.912856941806
```

The script will return the haplogroup and its score - the higher the better in terms of the quality of the assignment.
