With output files from angsd, bcf files were converted to .vcf. To merge all the .vcf files (one per chromosome), I had to use bgzip and tabix to merge the vcf files together

To run bgzip and tabix within my file directory:

```
for FILE in *.vcf; do bgzip $FILE; done

#source conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tabix

for FILE in *.vcf.gz; do tabix $FILE; done
```

After producing .vcf.gz files and index files, I then merged vcf.gz with bcftools
```
vcf-merge *.vcf.gz >  Salmo.salar.CTmax.combined.vcf
```
