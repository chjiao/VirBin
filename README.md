# VirBin   
Binning tools for strain-level viral contigs. The input to this program should be the contigs obtained from an assembly tool. We recommend to use PEHaplo (https://www.ncbi.nlm.nih.gov/pubmed/29617936) because it is opimized for assembling viral strains. Our tests also showed that PEHaplo's output can lead to better performance for viral contig binning. 

The output of this program are several sets of contigs. Each set corresponds to one haplotype (strain). The number of sets represent the number of haplotypes in the input data set. 

The method is essentially an EM algorithm. However, we created "windows" rather than whole contigs as input to the EM algorithm. The contigs have highly heterogeneous coverage because some parts are unique to a haplotype and some parts are common regions to multiple haplotypes. Our method employs sequence alignment between contigs to distinguish these cases. 

# Quick start
To quickly testing the core clustering algorithm, we have prepared the processed data set in folder *data/*.

## Dependencies
1. Install Python 2.7.x
2. Install Python module: networkx, numpy

## Installation 
Directly download the repository from github:   
```
git clone https://github.com/chjiao/VirBin.git
```

## Running example
You can test your program by running the following commands:
```
cd VirBin   
python VirBin.py -contig data/hiv5_contigs.fa -align data/hiv5_contigs_align.blastn -vcf data/hiv5_contig_align_0.9.vcf  -ref data/hiv5_ref_align.blastn   
```
If everything is good, you will see this information on your terminal
```
number of cluter: 5
length of the windows: 143
length of the groups duplicate windows: 68
Iteration: 4  
```
After running the commands above, there are four files output from the VirBin: contigs_windows.txt, windows_all.txt, preprocessed.blastn, EM_clusters.txt. The information of each windows is stored in contigs_windows.txt and windows_all.txt. The information of the final clusters is stored in EM_clusters.txt.

# The whole pipeline   
## Dependencies
1. Install Python 2.7.x
2. Install Python module: networkx, numpy
3. Install [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)    
4. Install [samtools](http://samtools.sourceforge.net/)   
5. Install [bcftools](https://samtools.github.io/bcftools/)   
Reminding: bowtie2, samtools and bcftools are recommended to be installed by [bioconda](https://bioconda.github.io/)   

## Pipeline steps
Input files: reads.fa, contig.fa, reference.fa   
1. Mapping reads on contigs   
Use script *tools/reads_align_on_contigs.sh*   
Output: Reads\_align\_on\_contig.vcf

2. Alignment between contigs
```
python ./tools/get_locations_single_blastn.py contig.fa contig_align.blastn
```

3. Align contigs on reference genomes (for evaluation only)   
```
python ./tools/get_locations_two_fas_blastn.py contig.fa reference.fa contig_ref_align.blastn
```
4. Run get_cluster_number.py (option)
If you **only** want to know the number of cluster, you can run the following command.
```
python VirBin.py -contig contig.fa -align contig_align.blastn -vcf Reads_align_on_contig.vcf -ref contig_ref_align.blastn
```

5. Run VirBin.py
The program will automatically decide the number of cluster and run the EM algorithm. Run as the following command.
```
python VirBin.py -contig contig.fa -align contig_align.blastn -vcf Reads_align_on_contig.vcf -ref contig_ref_align.blastn
```



