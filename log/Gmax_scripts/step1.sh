PREFIX=$1
THREAD=20
time bowtie2 -p "${THREAD}" -a -t \
	--end-to-end -D 20 -R 3 \
	-N 0 -L 10 -i S,1,0.50 --np 0 \
	--xeq -x index/gmax_rrna \
	-1 data/"${PREFIX}"/R1.fq.gz -2 data/"${PREFIX}"/R2.fq.gz | pigz -p 4 >output/"${PREFIX}"/rrna.raw.sam.gz
