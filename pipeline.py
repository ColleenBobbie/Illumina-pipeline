#Use PANDAseq to assemble paired-end Illumina reads into sequences. 
#Make directories to place the merged files in, and then use a 'for loop' to assemble the sequences

    mkdir /work/cbobbie/Illumina

    mkdir /work/cbobbie/Illumina/pandaseqmerged

    mkdir /work/cbobbie/Illumina/pandaseqmerged/outfiles
   

    for file in $(cat /work/cbobbie/samplenamescor.txt)
    do
      sqsub -q threaded -r 10m -n 4 -j Bacpandaseqmerge_${file} --mpp 2G -o /work/cbobbie/Illumina/pandaseqmerged/outfiles/${file}16Ssubmissionoutput.log 
      pandaseq -D 0.01 -f /work/cbobbie/allffiles/${file}_L001_R1_001.fastq.gz -r /work/cbobbie/allffiles/${file}_L001_R2_001.fastq.gz  -k 2 -l 390 -L 590 -p CCTACGGGAGGCAGCAG -q GGACTACHVGGGTWTCTAAT -t 0.6 -w /work/cbobbie/Illumina/pandaseqmerged/${file}.fasta -g /work/cbobbie/Illumina/pandaseqmerged/outfiles/${file}.log
    done

#combine all sequences into one fasta

    add_qiime_labels.py -i /work/cbobbie/Illumina/pandaseqmerged/ -m /work/cbobbie/map.txt -c InputFileName -n 1 -o /work/cbobbie/Illumina/combinedPANDA/

#make sure directory is in place for UPARSE commands

    mkdir /work/cbobbie/Illumina/uparse

#Dereplicate fasta
#usearch has a limit on sequences, so used derep_seqs.py from Gene Blanchard: https://gist.github.com/gblanchard4/8870470

     cd /work/cbobbie/Illumina/combinedPANDA 
     /home/cbobbie/panda-qiime/bin/packages/pipelinescripts/derep_seqs.py -i /combinedseqs.fna

     u8= /home/cbobbie/panda-qiime/bin/usearch8

#Abundance sort and discard singletons
    
    u8 -sortbysize /work/cbobbie/Illumina/uparse/derep_combined_seqs.fna -fastaout /work/cbobbie/Illumina/uparse/sorted.fasta -minsize 2

#OTU clustering

    u8 -cluster_otus /work/cbobbie/Illumina/uparse/sorted.fasta -otus /work/cbobbie/Illumina/uparse/otus1.fasta -relabel OTU_ -uparseout /work/cbobbie/Illumina/uparse/clustresults.txt
      "#4068 OTUs, 127667 chimeras (11.7%)"

#Chimera filtering
    usearch8 -uchime_ref /work/cbobbie/Illumina/uparse/otus1.fasta -db /work/cbobbie/SHARCNETsetup/work/databases/gold.fasta -strand plus -nonchimeras /work/cbobbie/Illumina/uparse/seqs.fna
      #Writing 3562 non-chimeras

#Map OTUS
    usearch8 -usearch_global /work/cbobbie/Illumina/combinedPANDA/combinedseqs.fna -db /work/cbobbie/Illumina/uparse/seqs.fna -strand plus -id 0.97 -uc /work/cbobbie/Illumina/uparse/otumap.txt
      #fatal error: File size too big for 32-bt version
      #8.5Gb
      
#Split original file into small chunks, then map OTUS.
    #First, create directories for the split and final files
        mkdir /work/cbobbie/Illumina/uparse/split_files
        mkdir /work/cbobbie/Illumina/uparse/uc_files
    #Then, download FASTA file splitter from http://kirill-kryukov.com/study/tools/fasta-splitter/
      #downloaded version 0.2.2 into /home/cbobbie/panda-qiime/bin/packages
    #In split_files, break the FASTA files into 100 equal parts:
      cd /work/cbobbie/Illumina/uparse/split_files
      perl /home/cbobbie/panda-qiime/bin/packages/fasta-splitter.pl  -n-parts-total 100 /work/cbobbie/Illumina/combinedPANDA/combinedseqs.fna
        #17333114 sequences dividing into 100 parts...
        "#All done, 451 seconds elapsed"
    #You will get the following files:
      ls /work/cbobbie/Illumina/uparse/split_files
    #Now run the following command in split_files folder: 
      for i in $(ls *.fna); do usearch8 -usearch_global $i -db /work/cbobbie/Illumina/uparse/seqs.fna -strand plus -id 0.97 -uc /work/cbobbie/Illumina/uparse/uc_files/$i.map.uc; done
        #This will map OTUs back to each split fasta
    #Check if uc files are there...
      ls /work/cbobbie/Illumina/uparse/uc_files
    #Finaly, combine all .uc files back together
      cd /work/cbobbie/Illumina/
      cat /work/cbobbie/Illumina/uparse/uc_files/* > ./map.uc

#create OTU table
    python /home/cbobbie/panda-qiime/bin/packages/pipelinescripts/uc2otutab_demultiplexed.py ./map.uc > ./otutable.txt

#assign taxonomy
    assign_taxonomy.py -i /work/cbobbie/Illumina/uparse/seqs.fna -r /work/cbobbie/SHARCNETsetup/work/databases/gg_97_otus.fasta -t /work/cbobbie/SHARCNETsetup/work/databases/gg_97_otu_taxonomy.txt -o /work/cbobbie/Illumina/assign_taxonomy/ -m rdp --rdp_max_memory 25000 -c 0.97

#convert to biom json formatted for use in phyloseq
    biom convert -i /work/cbobbie/Illumina/otutable.txt -o /work/cbobbie/Illumina/otutablenotax.biom --table-type ""OTU table"" --to-json

#add taxonomy info to biom
    biom add-metadata --sc-separated taxonomy --observation-header otu_table,taxonomy --observation-metadata-fp /work/cbobbie/Illumina/assign_taxonomy/seqs_tax_assignments.txt -i /work/cbobbie/Illumina/otutablenotax.biom -o /work/cbobbie/Illumina/illuminfinal.biom --output-as-json"

#convert Biom to txt (for easy access to admire all this work!)
   biom convert -i /work/cbobbie/Illumina/illuminfinal.biom -o /work/cbobbie/Illumina/illuminafinal.txt --to-tsv --header-key taxonomy

#Align sequences in QIIME, using greengenes reference sequences
#To build a real phylogenetic tree, we'll need to first align our sequences. We'll be using a reference alignment to align our sequences.
    align_seqs.py -i /work/cbobbie/Illumina/uparse/seqs.fna -o /work/cbobbie/Illumina/pynast_aligned -t /work/cbobbie/SHARCNETsetup/work/databases/core_set_aligned.fasta.imputed

#Filter alignments in QIIME
#This alignment contains lots of gaps from align_seqs.py, and it includes hypervariable regions that make it difficult to build an accurate tree. So, we'll filter it.  Filtering an alignment of 16S rRNA gene sequences can involve a Lane mask."
    filter_alignment.py -i /work/cbobbie/Illumina/pynast_aligned/seqs_aligned.fasta -o /work/cbobbie/Illumina/filtered_alignment -m /work/cbobbie/SHARCNETsetup/work/databases/lanemask_in_1s_and_0s

#Make the reference tree in QIIME
#Okay, let's make a tree out of that alignment!  This is actually quite easy in QIIME, using the make_phylogeny.py script (which uses the FastTree approximately maximum likelihood program, a good model of evolution for 16S rRNA gene sequences). The input for this script is our filtered alignment. We'll also root our tree at the midpoint, to ensure we can have a rooted tree for downstream analyses (i.e. Weighted UniFrac PCoA in phyloseq in R).
    make_phylogeny.py -i /work/cbobbie/Illumina/filtered_alignment/seqs_aligned_pfiltered.fasta -o /work/cbobbie/Illumina/rep_set.tre -r midpoint -t fasttree
