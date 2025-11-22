#!/bin/bash
#python gcf2omim.py --gcf GCFs/GCF_000001405.26_GRCh38_genomic.gff --mim mim2gene_clinical.txt > grch38_omim_clinical.bed
#python gcf2omim.py --gcf GCFs/GCF_000001405.25_GRCh37.p13_genomic.gff --mim mim2gene_clinical.txt > grch37_omim_clinical.bed

#run grch37
echo -e "chrom\tstart\tstop\tentrezGeneID\tGeneSymbol\tMIM" > ./grch37_omim_all_clinical.bed
python gcf2omim.py --gcf GCFs/GCF_000001405.25_GRCh37.p13_genomic.gff --mim mim2gene_clinical.txt | awk 'BEGIN { FS=OFS="\t" } {
    if ($2 == "") next
 
    # 1. Split the last column by "|" to handle multiple coordinates
    split($NF, coord_array, "|")
    
    # 2. Loop through each coordinate string found (will be 1 if no "|")
    for (j=1; j <= length(coord_array); j++) {
        
        # 3. Apply your original logic to the current coordinate string
        #    Split "chrY:59330367-59343488" into parts
        split(coord_array[j], parts, /[:\-]/)
        
        # 4. Print the new columns (chr, start, stop)
        printf "%s\t%s\t%s", parts[1], parts[2], parts[3]
        
        # 5. Print all original columns EXCEPT the last one (the duplication)
        for (i=1; i < NF; i++) {
            printf "\t%s", $i
        }
        
        # 6. Print the newline to finish this record
        printf "\n"
    }
}' >> ./grch37_omim_all_clinical.bed
cut -f1,2,3,5 ./grch37_omim_all_clinical.bed | sort | uniq > grch37_omim_clinical.bed

#run grch38
echo -e "chrom\tstart\tstop\tentrezGeneID\tGeneSymbol\tMIM" > ./grch38_omim_all_clinical.bed
python gcf2omim.py --gcf GCFs/GCF_000001405.26_GRCh38_genomic.gff --mim mim2gene_clinical.txt | awk 'BEGIN { FS=OFS="\t" } {
    if ($2 == "") next
 
    # 1. Split the last column by "|" to handle multiple coordinates
    split($NF, coord_array, "|")
    
    # 2. Loop through each coordinate string found (will be 1 if no "|")
    for (j=1; j <= length(coord_array); j++) {
        
        # 3. Apply your original logic to the current coordinate string
        #    Split "chrY:59330367-59343488" into parts
        split(coord_array[j], parts, /[:\-]/)
        
        # 4. Print the new columns (chr, start, stop)
        printf "%s\t%s\t%s", parts[1], parts[2], parts[3]
        
        # 5. Print all original columns EXCEPT the last one (the duplication)
        for (i=1; i < NF; i++) {
            printf "\t%s", $i
        }
        
        # 6. Print the newline to finish this record
        printf "\n"
    }
}' >> ./grch38_omim_all_clinical.bed
cut -f1,2,3,5 ./grch38_omim_all_clinical.bed | sort | uniq > grch38_omim_clinical.bed
