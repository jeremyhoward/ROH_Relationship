# ROH_Relationship

The ROH_Relationship is a program that creates a Run of Homozygosity (ROH) based relationship matrix. A full description of how it is calculated is outlined in "Characterization and management of long runs of homozygosity in parental nucleus lines and their associated crossbred progeny".

## To Compile and Run
Compile: g++ -O3 -fopenmp -lpthread ROH_Rel.cpp -o ROH_Relationship

Run: ./ROH_Relationship Test_Genotype Test_Map 5 2 4 ROH5

## Parameter File
The program takes in 6 required parameters and some information about them is outlined below:

- **Genotypes:** The name of the genotype file and has to be a single word. The format is ID and followed by genotype string. The string of genotypes should be the same length across individuals. The delimiter is a space. The genotypes can’t be missing and they have to be in the phased format 0 (homozygote (AA)), 3 (heterozygote Aa), 4 (heterozygote aA), and 2 (other homozygote (aa)). The row of the map should correspond the the location within the genotype string.
- **Map File:** The name of the map file and has to be a single word. The format of the map file is chromosome for column 1 and nucleotide position for column 2. The map file has to be ordered by chromosome and nucleotide position within chromosome. The delimiter is a space.
- **ROH Cutoff:** Refers to how long your ROH stretches have to be in nucleotides in order to be called a ROH stretch. For example a value of “5000000” refers to a cutoff of 5 Mb.
- **ROH Threshold:** This threshold is given to remove regions that have a small number of SNP and therefore could give rise to spurious results. The value refers to how many standard deviations from the mean number of SNP that are removed. For instance a value of 2 refers to any window that contains less than 2 * Average SNP will be removed.
- **Threads:** Number of cores to utilize when constructing relationship matrix.
- **Out File:** Name of file to output relationship. Has to be one word. Relationships are outputted in triplet form. 
