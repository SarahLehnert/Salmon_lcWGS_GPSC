This script is for after intial pre-processing and genotyping. This script includes filtering the data and running imputation/phasing, as well as cleaning up individual ID names, adding SNP names, updating chromosomes numbers, adding phenotype and population info.


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

Also need to update Chromosome numbers (currently names from RefSeq) - so update as numbers instead, and also update SNP names (currently no SNP names assigned). See chromosome numbers at https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Salmo_salar/102/

To update these, I ran an R script to replace Chr numbers and create SNP ID then overwrite the .bim plink file with the updated info.

```
#Update Chromosome # and assign SNP names

setwd("/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/")

Chr_names<-read.table("Chromosome_Name_Update.txt", header=T)
bim_file<-read.table("Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.bim", header=F)

head(Chr_names)
head(bim_file)

#create column with order to keep the same order on output
bim_file$Order<-c(1:nrow(bim_file))

combined<-merge(x=bim_file, y=Chr_names, by.x="V1", by.y="RefSeq")

#check that same number of SNPs remain
nrow(combined)
nrow(bim_file)

#update SNP names (column V2) to be "Chr#_Position"
combined$V2 <- interaction(combined$Chr, combined$V4, sep="_")

#check
head(combined)

outputFile <- combined[order(combined$Order),]
head(outputFile[,c(8,2:6)])


#now make sure order is correct, and save to overwrite original .bim file - this should work to update Chr and SNP names
write.table(outputFile[,c(8,2:6)], "Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.bim",quote = F, row.names = F, col.names = F, sep="\t" )

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

#Create new sample ID name
newd$SampleName <- paste0("CTmax_", newd$V1)

#Update any unknown pops as "unknown" instead of different levels here:
newd$Type.1 <- as.character(newd$Type.1)
newd$Type.1[which(newd$Type.1=="?")] <- "unknown"
newd$Type.1[which(newd$Type.1=="blu yell?")] <- "unknown"
newd$Type.1[which(newd$Type.1=="CMP?")] <- "unknown"

#Re-order to original order
newd<-newd[order(newd$order),]
head(newd)

#open original fam file to get original names
fam_org<-read.table("/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.fam", header=F)
head(fam_org)


##Could create a new file to run with plink and update ID names... or could just overwrite .fam file here:

#Create file needed to update-ids function in plink
out<-as.data.frame(cbind(as.character(fam_org$V1),
      as.character(fam_org$V2),
      as.character(newd$SampleName), as.character(newd$Type.1)))

#save file and run plink to update names
write.table(out, file = "/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/update_IDs_ctmax.txt", quote = F, row.names = F, col.names = F, sep="\t")



#Save file to overwrite fam:


out_fam<-as.data.frame(cbind(as.character(newd$SampleName), 
                             as.character(newd$Type.1), 
                             as.numeric(newd$V3),as.numeric(newd$V4),as.numeric(newd$V5),
                             as.numeric(newd$Temp) ))

#Note - this will overwrite file - so won't be able to re-run this script again - unless running previous step in plink again...
write.table(out_fam, file = "/Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final.fam",
            quote = F, row.names = F, col.names = F, sep="\t")


````

Next use output file with updated info to re-run plink 
Also run plink with --recode 12 and --recode A for files that may be needed for GWAS

```
./plink --bfile /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final \
 --allow-extra-chr --chr-set 29 \
 --make-bed --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final_updateID

./plink --bfile /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final \
 --allow-extra-chr --chr-set 29 \
 --recode 12 --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final_updateID

./plink --bfile /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final \
 --allow-extra-chr --chr-set 29 \
 --recode A --out /Volumes/Accelsior_8TB/Salmon/CTmax_genotype/Salmo.salar.CTmax.combined_maf005_geno07_phased_imputed_Final_updateID


```

##also ran plink with --recode 12 and --recode A for files that may be needed for GWAS
#next trying pcadapt
