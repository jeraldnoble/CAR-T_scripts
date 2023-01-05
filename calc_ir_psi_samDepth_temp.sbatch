#!/bin/bash
#SBATCH --job-name="BASE_irPSI"
#SBATCH --output=BASE_irPSI.out
#SBATCH --error=BASE_irPSI.err
#SBATCH --mail-type=FAIL
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --time=200:00:00
#SBATCH --no-requeue
#SBATCH --mail-user="jerald.noble@moffitt.org"

cd $SLURM_SUBMIT_DIR

module load BCFtools
module load SAMtools
ml python3

bed=/share/lab_robinson/jd/carT/patient_cohort/CD19_gene_coords.bed

samtools view -h BASE_filt.bam -L ${bed} > BASE_cd19.sam

python3 /share/lab_robinson/jd/scripts/filt_chimeric_reads_from_CD19_reads.py -sam BASE_cd19.sam -bed ${bed} > BASE_chim_filt.sam

samtools depth BASE_chim_filt.sam > BASE_cd19_depth.tsv

python3 /share/lab_robinson/jd/scripts/calc_intron_psi_from_samDepth.py -depth BASE_cd19_depth.tsv -sj BASE_filt_sj.tsv \
	-bed /share/lab_robinson/jd/carT/patient_cohort/CD19_ir_coords.bed > BASE_cd19_ir.tsv

rm BASE_cd19.sam BASE_chim_filt.sam BASE_cd19_depth.tsv
