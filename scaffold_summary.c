#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "fasta.h"

int intcomp(const void *x,const void *y){
    return *(int *)y-*(int *)x;
}



int main(int argc,char **argv)
{
    FASTAFILE *ffp;
    char *seq;
    unsigned long fasta[10000];
    char *seqname,*L_scaffold_name;//name of current sequence
    int seqlen;   //length of current sequence
    unsigned long L_scaffold_num = 0;   //numbers of lagrge scarffold (1000bp+)
    unsigned long L_scaffold_base = 0;  //bases of large scarffold (1000bp+)
    unsigned long total_base = 0;    //total bases of all scaffold
    unsigned long L_scaffold_len = 0;   //length of largest scaffold
    unsigned long n50 = 0;
    unsigned long n90 = 0;
    unsigned long sum = 0;
    int n50_num,n90_num;
    int i,n = 0;
    unsigned long total_seq = 0;
    float num_n1 =0;
    float num_n2 =0;
    float num_gc1 =0;
    float num_gc2 =0;
    float gc_content1,n_ratio1;
    float gc_content2,n_ratio2;
	
    if(argc==1||argc>2){
	printf("\nUsage: scaffold_summary <*.fa>\n\tBrief summary of N50,N90 etc. of sequences in specificed fasta.\n\tWenchao Lin <liwnenchao@yeah.net>\n\n");
	exit(0);
    }

    ffp=OpenFASTA(argv[1]);
    while(ReadFASTA(ffp,&seq,&seqname,&seqlen))
    {
        fasta[n]=seqlen; 
        n++;
        if(seqlen>L_scaffold_len){
            L_scaffold_len = seqlen;
            L_scaffold_name = seqname;
        }
        if(seqlen>1000){
            L_scaffold_num ++;
            L_scaffold_base += seqlen;
            for (i=0;i<seqlen;i++)
            {
                if(seq[i]=='g' || seq[i]=='G'||seq[i]=='c'||seq[i]=='C')
                {
                    num_gc1++;
                }
                if(seq[i]=='n' || seq[i]=='N')
                {
                    num_n1++;
                }
            }
        }
        for (i=0;i<seqlen;i++)
        {
            if(seq[i]=='g' || seq[i]=='G'||seq[i]=='c'||seq[i]=='C')
            {
                num_gc2++;
            }
            if(seq[i]=='n' || seq[i]=='N')
            {
                num_n2++;
            }
        }

        total_base += seqlen;

        free(seq);
    }
    qsort(fasta,n,sizeof(unsigned long),intcomp);
    //    printf("Ns: %1.0f\n",num_n);
    //    printf("GCs: %1.0f\n",num_gc);
    //    printf("total: %d\n",n);
    total_seq = n;
    gc_content1 = num_gc1*100/L_scaffold_base;
    n_ratio1 = num_n1*100/L_scaffold_base;
    gc_content2 = num_gc2*100/total_base;
    n_ratio2 = num_n2*100/total_base;

    for(i=0;i<n;i++){
        sum += fasta[i];
        if(sum>=0.5*total_base && n50 == 0){n50 = fasta[i];n50_num = i+1;}
        if(sum>=0.9*total_base && n90 == 0){n90 = fasta[i];n90_num = i+1;break;}
    }
    printf("File: %s\n",argv[1]);
    printf("Large scaffolds (>1000bp)\n");
    printf("\tLargest scaffold: %s\n\tLargest length: %ld\n",L_scaffold_name,L_scaffold_len);
    printf("\tScaffolds 1000bp+: %ld\n\tBases in largest scaffold: %ld\n",L_scaffold_num,L_scaffold_base);
    printf("\tN50 scaffold: %d\n\tN90 length: %ld\n",n50_num,n50);
    printf("\tN90 scaffold: %d\n\tN90 length: %ld\n",n90_num,n90);
    printf("\tG+C content: %f%%\n\tN rate: %f%%\n",gc_content1,n_ratio1);
    printf("\nAll scaffolds\n\tNumber: %ld\n\tTotal bases: %ld\n\tGC content: %f\n\tN rate: %f\n",total_seq,total_base,gc_content2,n_ratio2);
   CloseFASTA(ffp);
    exit(0);
}
