module load SAMtools
ml python3

bed=CD19_gene_coords.bed

### start here with your BAM file
### BASE is the sample name, needs to be replaced for a given sample via sed or tr

samtools view -q 30 -F 256 --threads 3  -bh BAM  > BASE_filt.bam

samtools view BASE_filt.bam | awk -f /share/lab_robinson/jd/scripts/sjFromSAMcollapseUandM.awk | sort -V > BASE_filt_sj.tsv

samtools view -h BASE_filt.bam -L ${bed} > BASE_cd19.sam

python3 filt_chimeric_reads_from_CD19_reads.py -sam BASE_cd19.sam -bed ${bed} > BASE_chim_filt.sam

samtools depth BASE_chim_filt.sam > BASE_cd19_depth.tsv

python3 calc_intron_psi_from_samDepth.py -depth BASE_cd19_depth.tsv -sj BASE_filt_sj.tsv \
	-bed CD19_ir_coords.bed > BASE_cd19_ir.tsv

rm BASE_cd19.sam BASE_chim_filt.sam BASE_cd19_depth.tsv
