### The goal of this Unix shell script is to filter the SNP data by MAP using VCFTools and get them ready for sNMF analyses.

# Filtering SNPs by MAF and then extracting one SNP per locus for SNMF analyses
# From Ipyrad vcf format output file

# Select our data (name based on ipyrad paramaters)
for s in fusco_R12_c95_n162_m120 #fusco_R12_c95_n158_m120
do

# create a folder and change directory to that new folder
mkdir /Users/dequeirk/Dropbox/2018_fuscoauratus/VCFtools_SNMF
mkdir /Users/dequeirk/Dropbox/2018_fuscoauratus/VCFtools_SNMF/${s} 
cd /Users/dequeirk/Dropbox/2018_fuscoauratus/VCFtools_SNMF/${s}

# Filter by MAF = 0.05 using VCF Tools
vcftools --vcf /Users/dequeirk/Dropbox/2018_fuscoauratus/ipyrad_outfiles/${s}_outfiles/${s}.vcf --recode --out ${s}_filtered --min-alleles 2 --max-alleles 2 --maf 0.05

# Now extracting one random SNP per locus
cat ${s}_filtered.recode.vcf | grep "#" > ${s}_1SNP-locus.vcf # extract the headers (whose line starts with #), save in different file
cat ${s}_filtered.recode.vcf | grep -v "#" | sort -u -k1,1 >> ${s}_1SNP-locus.vcf # sort by 'chromosome' name and extract unique, concatenate to that file created

# Now saving in 012 format
vcftools --vcf ${s}_1SNP-locus.vcf --out ${s}_1SNP-locus --012

# Now replace -1 with 9 to indicate missing sites
cat ${s}_1SNP-locus.012 | sed -e 's/-1/9/g' ${s}_1SNP-locus.012 > ${s}_1SNP-locus_tmp.012 # replace and save in temporary file
mv ${s}_1SNP-locus_tmp.012 ${s}_1SNP-locus.012 # rename temp file as final 012 file

# close the loop:
done

# Done!

