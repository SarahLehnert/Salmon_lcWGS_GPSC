
First filter out individuals with low genotyping rate (<30%) - current file includes some blanks so these should be removed before imputation/phasing. 
Here we used plink to filter out individuals with <30% genotyping rate, and also filtered for MAF>0.05

```
cd /Desktop/Software/plink_mac_20200219
./plink --vcf /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined.vcf.gz --allow-extra-chr \
--double-id --mind 0.3 --maf 0.05 \
--recode vcf --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno03.vcf
```

Using the filtered dataset, we then ran Beagle to phase/impute data

```
cd /Desktop/Software
java -Xmx300G -jar beagle.18May20.d20.jar gt=/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno03.vcf.vcf out=/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed  nthreads=30
```

Phased/imputed genotype data (vcf) was then converted into plink files for further analyses:

```
cd /Desktop/Software/plink_mac_20200219

./plink --vcf /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed.vcf.gz \
 --allow-extra-chr --maf 0.05 --double-id \
 --make-bed --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final
```

Now IDs need to be updated - remove duplicated part within IDs and additional part related to sequencing info. Also added "Population" identifier (WAB, CMP, GAR) for individuals. 

Updating IDs was done through R using a script to create updated ID file to run in plink.

```
#read in fam file
fam<-read.table("/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.fam", header=F)

#remove extra character strings in column 1 names
fam$V1<-gsub(pattern = ".*CTmax", replacement="", x =fam$V1 )
fam$V1<-gsub(pattern = ".*CTMax", replacement="", x =fam$V1 )
fam$V1<-gsub(pattern = ".realigned.bam", replacement="", x =fam$V1 )
fam$V1  <-as.numeric(as.character(fam$V1))

#remove extra character string in column 2 names
fam$V2<-gsub(pattern = ".*CTmax", replacement="", x =fam$V2 )
fam$V2<-gsub(pattern = ".*CTMax", replacement="", x =fam$V2 )
fam$V2<-gsub(pattern = ".realigned.bam", replacement="", x =fam$V2 )
fam$V2  <-as.numeric(as.character(fam$V2))

#create a column with "order" to keep in same order
fam$order <- as.numeric(c(1:nrow(fam)))

head(fam)

#read in ctmax data - with info for IDs
ctmaxdat<-read.csv("/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/CTmax_Data_trials.csv", header=T)

#merge fam and ID info
newd<-merge(fam, ctmaxdat[,c(1,6,10)], by=1)

#Re-order to original order
newd<-newd[order(newd$order),]

#open original fam file to get original names
fam_org<-read.table("/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.fam", header=F)

#Create file needed to update-ids function in plink
out<-as.data.frame(cbind(as.character(fam_org$V1),
      as.character(fam_org$V2),
      as.character(newd$V1), as.character(newd$Type.1)))

#save file and run plink to update names
write.table(out, file = "/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/update_IDs_ctmax.txt", quote = F, row.names = F, col.names = F, sep="\t")
````

Next use output file of updated IDs to run with plink 

```
./plink --bfile /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final \
 --allow-extra-chr --update-ids /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/update_IDs_ctmax.txt \
 --make-bed --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final_updateID
```

##also ran plink with --recode 12 and --recode A for files that may be needed for GWAS
#next trying pcadapt
