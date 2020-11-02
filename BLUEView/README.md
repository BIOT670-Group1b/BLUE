## OVERVIEW:

* RunBLUE.sh is the bash script that runs the application
    * To run: Open linux terminal, set directory to checkout folder and execute command ./RunBLUE.sh
        * You may encounter missing package. Install any that are not found
* BLUE comes with sample BAM files to test with in the "Sample BAM" folder
* Application folder contains the following:
    * Python Applciation Scripts
    * Chromosomes folder
        * Place for the Homo_sapiens.GRCh38.dna.chromosome fils
        * They are too big and need to be download and place in this folder manually
            * Download
            * wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz
            * wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{X,Y}.fa.gz
	* Genomes (GTF Reference file(s))
        * Place for the Homo_sapiens.GRCh38.100.db and Homo_sapiens.GRCh38.100.gtf
        * They are too big and need to be download and place in this folder manually
            *  Download and unzip
            * wget ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
            * gunzip Homo_sapiens.GRCh38.100.gtf.gz
	
## USING THE APPLICATION:

### BAM ONLY: 
To use this application, simply select a BAM file from the 'File' dropdown and the program will read it. If a fusion is found, the reads will be displayed in the Fusion Area in the middle of the application. The first fusion gene will be displayed and color coded on the left side, while the second will show up on the right. Positions mark the start of the base that follows and each 'marker' is 5 bases in length.

### BAM W/ ANNOTATIONS: 
In order for the program to retrieve the gene names and other relevant information, a GTF reference file will need to be loaded into the program. This can be done by selecting the reference the 'Ref' tab and the reference name associated with that genome. The program comes preloaded with Ensembles GChr38 annotation file for use. You can load the reference file before loading the BAM file and vice versa. For fusions found, the gene list A and the gene list B will populate with discovered fusions genes. More information about each gene can be found by clicking 'More Information' under each gene name. The summary section will show key information about the fusion that was found.

### USING THE CHROM-VIEW:
The program allows each chromosome sequence to be viewed on a viewer under the Fusion-View. To load a chromosome sequence, select any of the chromosomes from the dropdown near the top of the window. The area below the Chromosome viewer updates in real-time the genes currently being viewed in the viewer and their positions. For the genes to be updated, the reference file must be loaded into the program, otherwise no genes will populate. Genes and positions can also be searched by using the search function next the the chromosome dropdown list. Type in the gene name or single position and the chrom-view will automatically go to the specified region/gene if it exists. The Fusion-View will update the same way, though will not show any sequences outside of the fusion sequences.
