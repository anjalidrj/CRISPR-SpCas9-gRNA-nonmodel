# If your .gff3 does not have the target gene annotated
## Steps to create a bed file on notepad:
1. Find the position of the gene in the genome sequences. Specifically, find start and end points. This can be done using snapgene, benchling, etc.   
2. Open a blank txt file in notepad.  
3. Type in the chromosome/sequence name (will be given in the fasta file, or can be found when generating other bed files for the organism).  
4. Click on tab 
5. Then type in starting nucleotide number.
6. Click on tab then type in ending nucleotide number.
7. Click on tab and type in your gene's name.
8. Click on tab and type 0.
9. Click on tab and type either + or - depending on whether the sequence is on the top or bottom strand.
10. At the end, the file should look like this:
    U00096.3	10642	11356	yaaW	0	-
11. Save it as gene_name.bed 
