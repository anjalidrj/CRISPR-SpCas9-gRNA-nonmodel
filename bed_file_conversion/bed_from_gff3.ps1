#Run one at a time in powershell(windows)

#1
cd "C:\Users\name\folder_loaction\folder_name" # <--- copy folder containing fasta and gff3 files pathname

#2
# Preview first 3 lines matching "gene_name"
# Change annotation.gff3 to your .gff3 file's name and "gene_name" to the name of your gene
Get-Content annotation.gff3 | Select-String "gene_name" | Select-Object -First 3

#3
# Extract gene entry and convert to BED format
# Change "gene_name" in the 4th, 6th and 8th line and change annptation.gff3 in the 1st line to your .gff3 file's name.
Get-Content annotation.gff3 |  
    Where-Object { $_ -notmatch "^#" } |
    ConvertFrom-Csv -Delimiter "`t" -Header seqid,source,type,start,end,score,strand,phase,attributes |
    Where-Object { $_.type -eq "gene" -and $_.attributes -match "Name=gene_name" } |
    ForEach-Object {
        "$($_.seqid)`t$([int]$_.start - 1)`t$($_.end)`tgene_name`t0`t$($_.strand)"
    } |
    Out-File -Encoding ascii "gene_name.bed"

#4
# Verify output, change the gene_name to match your gene's name. 
Get-Content gene_name.bed
