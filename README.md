 gcf2omim
Convert NCBI RefSeq/GenBank GCFs to OMIM disease-related gene coordinates 
Based on the hg19Generator script provided by OMIM, modified to convert any GCF to OMIM using mim2gene.txt

Required files:
GRCh37 1405.25
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz

GRCh38 1405.26
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.gff.gz

And mim2gene.txt which can be downloaded from OMIM:
https://omim.org/downloads

Example scripts: 
Generating disease-associated OMIM gene coordinates for both hg19 and hg38: 
bash get_omim_clinical.sh
Mapping gene overlaps between a given set of genomic coordinates to hg38 OMIM genes: 
bash example_intersect.sh my_genomic_regions.bed grch38_omim_clinical.bed
#Summarizing two sets of genomic coordinates intersected with hg19 and hg38 separately into a TSV comparing concordance in genes annotated 
python merge_counts.py --hg19 my_hg19_omim_overlaps.bed --hg38 my_hg38_omim_overlaps.bed
