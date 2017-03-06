# annotation-scripts
Useful scripts for utilizing BLAST/Exonerate and manipulating resulting output and GFF files.

These were employed in the annotation of the fathead minnow genome assemblies.

## Contents
File | Purpose
---- | ----
exonerate_afterblast.pl | Coordinates the use of BLAST+ and Exonerate to produce spliced alignments on a relatively fragmented genome assembly.
multithread_exonerate_afterblast.pl | Asynchronous version of exonerate_afterblast for speed improvement
\*.config files | Configuration files supplied to exonerate_afterblast programs, can be customized for different datasets
simple_exonerate2gff3.pl | easy method for converting exonerate output (with GFF output enabled) to GFF3 format
FilterGff.pl | Filter GFF3 entries using customized criteria
bestOfGFF.pl | Select best of GFF3 entries for a particular query sequence
add_tags.R | use BiomaRt to add gene symbols and descriptions to GFF3 files


## Usage
The entire process is outlined for the alignment of fathead minnow (*Pimephales promelas*) Expressed Sequence Tags (ESTs) to the SOAPdenovo genome assembly [(Burns et al. 2015; Environ. Toxicol. Chem. 35:212)][1]. This process can be extended to other datasets as desired.

### Retrieving sequences

The SOAPdenovo assembly was downloaded from NCBI's FTP site and decompressed:

    #Create directory for analysis and cd into it
    mkdir ~/FHM_Genome
    cd ~/FHM_Genome/
    curl -o SOAPdenovo.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Pimephales_promelas/all_assembly_versions/GCA_000700825.1_FHM_SOAPdenovo/GCA_000700825.1_FHM_SOAPdenovo_genomic.fna.gz
    gzip -d SOAPdenovo.fna.gz

A Multifasta of fathead minnow ESTs were downloaded from NCBI's nucest database using the web browser with the following address:

    https://www.ncbi.nlm.nih.gov/nucest/?term=pimephales+promelas+%5BOrganism%5D

The resulting file was copied into the previously created FHM_Genome directory.

Each of these multifasta files (genome and ESTs) were split into single-entry fasta files for use in the next steps. [A simple bioPerl program][2] was used for this, but many options are available for this task. The end result was to have two new directories - one containing 258,504 individual EST query fasta files, and another containing 73,058 individual FHM_SOAPdenovo scaffolds from the genome assembly.

Additionally, a blast database was created for the fathead minnow genome assembly:

    makeblastdb -dbtype 'nucl' -in SOAPdenovo.fna -title 'FHM_assembly' -parse_seqids -hash_index -out 'FHM_assembly' -taxid 90988

### Aligning sequences

The program exonerate_afterblast.pl was used to produce spliced alignments using BLAST+ and Exonerate. The EST configuration file 'EST.config' from this repository contains alignment settings for each program.

    mkdir ~/FHM_Genome/EST_aln_out/
    cd ~/FHM_Genome/EST_aln_out/
    #Make sure to set BLASTDB variable before running
    #e.g. export BLASTDB=/path/to/blast_database/
    perl exonerate_afterblast.pl -c ~/FHM_Genome/EST.config -q ~/FHM_Genome/EST_queries/ -t ~/FHM_Genome/Subject_scaffolds/ 2>&1 | tee ../alignment_log.txt

The result is a directory full of exonerate output files, one for each query sequence.

### Processing Exonerate output

Exonerate's output was configured (via the config file) to contain GFF formatted alignment information. However, this older version of GFF format contained in the output is deprecated and connot encode parent-child relationships (e.g. gene->mRNA->exon structure). This output can nevertheless be converted into GFF3 format (which can encode parent-child relationships), and the program simple_exonerate2gff3.pl was employed for this purpose. Note that in this example I concurrently process customized roll-your-own exonerate output in order to include alignment coverage information for each completed alignment. This information will be used in the next step.

    #Results were concatenated into single file before processing
    find ~/FHM_Genome/EST_aln_out/ -type f -iname "*.exonerate_out" | xargs cat > ESTs_all.exonerate_out
    #Converted the exonerate output to GFF3 formatted file
    perl ~/annotation-scripts/simple_exonerate2gff3.pl -P -R '^---END RYO---' ESTs_all.exonerate_out > ESTs_all.gff3

### Manipulating GFF3 entries

After the data had been converted to GFF3 format, the alignment entries which covered less than 90 percent of the query sequence were removed.

    perl ~/annotation-scripts/FilterGff.pl -f 'aln_pct=(.*?);' -v 90 ~/FHM_Genome/ESTs_all.gff3 > ~/FHM_Genome/ESTs_filtered_aln_pct_90.gff3

Next, all aligment entries except the highest-scoring for each query were removed.

    perl ~/annotation-scripts/bestOfGff.pl ~/FHM_Genome/ESTs_filtered_aln_pct_90.gff3 > ~/FHM_Genome/ESTs_best.gff3

The resulting GFF3 file contained many rows of redundant information, and so the result was simplified to contain only the informative 'gene' or 'exon' type entries

     awk '($3=="gene"||$3=="exon"||$0=="###")' ~/FHM_Genome/ESTs_best.gff3 > ESTs_best-min.gff3

Now we have created a gff3 file which describes the highest-scoring alignments between the EST query sequences and the fathead minnow SOAPdenovo genome assembly. This file can be viewed in a genome browser or it can be used to simply move this information into other bioinformatic analyses.

[1]: https://www.ncbi.nlm.nih.gov/pubmed/26513338
[2]: http://www.perlmonks.org/?node_id=767502
