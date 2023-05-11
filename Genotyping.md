We will calculate [genotype likelihoods](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3593722/), which are uncertainty-weighted estimates of an individual's genotype at a locus. We will carry this out across each chromosome separately, because it is somewhat memory and time intensive. Notice that we are now exporting a couple of different variables to slurm - we are specifying chromosome rather than set, and some run-specific parameters that we keep in another params file.

Here will need a list of chromosomes names in input.

```
while read chrom;  do sbatch --export=ALL,chrom=$chrom,paramfile=WGSparams_CTmax.tsv,angsdparam=refs_angsdparam.tsv  09_angsd_bcf_beag_maf.sh ;  
  done < ../Ssal_v3.1_genomic.chroms
```
For the CTmax project, I ran with all individuals together (not reference samples).

The ANGSD command looks like this:

```
cd projdir/angsd_in

$angsd \
  -nThreads 8 \
  -bam $bamfile \
  -out $projdir/angsd_out/$species.$projname.$runname.$chrom. \
  -dobcf 1 \
  -gl 1 \
  -dopost 1 \
  -dogeno 5 \
  -doGlf 2 \
  -domajorminor 1 \
  -domaf 1 \
  -docounts 1 \
  -dumpCounts 2 \
  -doQsDist 1 \
  -minMapQ 30 \
  -minQ 30 \
  -minInd $minInd \
  -SNP_pval 2e-6 \
  -uniqueOnly 1 \
  -minMaf 0.05 \
  -setMinDepth $minDepth \
  -r $chrom \
  -remove_bads 1 \
  -only_proper_pairs 1
 ```
Most of it won't change, but we can use sites genotyped at high coverage for better imputation (potentially? We'll see). To get these sites after the first genotyping run, we do the following:

```
zcat *mafs.gz  | sort | uniq | cut -f1,2,3,4 | sed '$d' > All_sites.tsv
```

Which will open up all the allele frequency estimates from each chromosome, sort them and remove duplicate headers, and then drop the header that gets sorted to the bottom of the file.  Then we index these sites using ANGSD:

```
conda activate align
angsd sites index All_sites.tsv
conda deactivate 
```
We can add this information to the parameter file for our lcWGS samples (as below) - however, I did only ran script 09, and not 10, as I did not use a set of reference sites/individuals. 

For minInd use 0.80 x # of individuals - indicating SNP is genotyped in >80% of samples
For minDepth use 2 x # of individuals - indicating ~2x coverage for genotyping

```
#Update parameter file 
bamfile=refbam_Ctmax.tsv
runname=refsamples
minInd=307
minDepth=614
sites=angsd_in/All_sites.tsv

```
Note I did not run 10_refsites_angsd_bcf_beag_maf.sh

Now we have a bunch of kinda useless bcf files. They're smaller than vcfs but they lack human interpretable info, so we need to convert them to vcfs. We do this per bcf file, so we again specify the chromosome and file set.

```
while read chrom;  do sbatch --export=ALL,chrom=$chrom,paramfile=WGSparams_aeip.tsv,angsdparam=refs_angsdparam.tsv  11_bcf_to_vcf.sh ;  
  done < Ssal_v3.1_genomic.chroms
```
Next I combined vcf files. To merge all the .vcf files (one per chromosome), I had to use bgzip and tabix to merge the vcf files together.

To run bgzip and tabix within my file directory:

```
for FILE in *.vcf; do bgzip $FILE; done

#source conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tabix

for FILE in *.vcf.gz; do tabix $FILE; done
```

After producing .vcf.gz files and index files, I then merged vcf.gz with bcftools concatenate

```
#vcftools seems to be taking a long time - trying bcftools instead

conda activate bcftools
bcftools merge *vcf.gz > Salmo.salar.CTmax.combined.vcf
deactivate bcftools

bgzip Salmo.salar.CTmax.combined.vcf

conda activate tabix
tabix Salmo.salar.CTmax.combined.vcf.gz

```
Save .vcf.gz to local computer
```

scp username@address:/filepath/Salmo.salar.CTmax.combined.vcf.gz IP_Address:C:/localdirectory/
