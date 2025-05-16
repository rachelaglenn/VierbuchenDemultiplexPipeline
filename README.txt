#1. Run the create genome script (1_create_genome.txt) to generate the mouse reference genome used by the dropseq pipeline and STAR aligner. You only have to run this step once and you can use the reference genome for all samples. 

#2. Concatenate all lanes for each sample using cat (see example script concat_fastqs.txt for sample C)

#3. Run the Dropseq pipeline script (dropseq_alignment_script.sh). The script takes 3 agruments: fastq1, fastq2 and sample name. The input path is specified in the script, you can change the line accordingly. 
bash 3_dropseq_alignment_script.sh C_P2-5-11-25-30_merged_R1_001.fastq.gz C_P2-5-11-25-30_merged_R2_001.fastq.gz C_P2-5-11-25-30_merged

#4. Run the detect doublets script (4_detect_doublets_all.txt) to calculate the donor assignment for each cell. 

#5. Run the R script (5_own_cell_assignment.R) for the updated donor assignment method.

#6. Once all samples were processed using step 2 - 5, run the integration script (6_integration.R) for downstream analysis. 
