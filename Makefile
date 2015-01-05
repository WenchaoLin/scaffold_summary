objects = scaffold_summary.o fasta.o
scaffold_summay : $(objects)
	gcc -o scaffold_summary $(objects)
scaffold_summary.o :  scaffold_summary.c fasta.h
	gcc -c scaffold_summary.c
fasta.o : fasta.c fasta.h
	gcc -c fasta.c
clean : 
	rm scaffold_summary $(objects)
