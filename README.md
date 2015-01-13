# Synopsis

A simple program for evaluate the assembly of a give genome. Shows simple statistics about a fasta file containing many sequences. The result includes the name, N50, N90, GC Content, N rate, number of sequences etc., as well as various other information both in large scaffolds (>1000 bp) and in all sequences. Command line usage for scaffold_summary can be found by typing scaffold_summary. 

# Usage

Command line usage can be found by running the following:

`scaffold_summary <*.fa>`

# Installing

Complie the source code by typing `make` in the  folder.

# Sample output


```
File: 454LargeContigs.fna

Large scaffolds (>1000bp)
	Largest scaffold: contig00001 length=24332   numreads=1033
	Largest scaffold size: 24332
	The number of large Scaffolds: 154620
	Bases in large scaffolds: 351138116
	The number of N50 scaffolds: 56746
	N50 length: 2253
	The number of N90 scaffolds: 159448
	N90 length: 960
	G+C content: 4.7780%
	N rate: 0.0000%

All scaffolds
	The number of sequences: 214338
	Total bases: 395414334
	GC content: 4.2429%
	N rate: 0.0000%
```

# Author

	Wenchao Lin <linwenchao@yeah.net>
	Tianjin Biochip Corporation @ 2015
