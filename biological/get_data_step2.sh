mkdir data;cd data;
mkdir geo;cd geo;

# https://doi.org/10.1038/s41586-022-05580-6
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM4_ESM.xlsx

# https://doi.org/10.1038/s41467-020-20225-w
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20225-w/MediaObjects/41467_2020_20225_MOESM3_ESM.xlsx

wget -O Human.GRCh38.p13.annot.tsv.gz https://www.ncbi.nlm.nih.gov/geo/download/?format=file\&type=rnaseq_counts\&file=Human.GRCh38.p13.annot.tsv.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210351/matrix/GSE210351_series_matrix.txt.gz

wget -O GSE210351_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts\&acc=GSE210351\&format=file\&file=GSE210351_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz

wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt
wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz
gzip -d Homo_sapiens.GRCh38.111.gtf.gz 

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222147/suppl/GSE222147_lcm_wbgs_bsseq_smoothed_coverage_filtered_geo.rda.gz;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142241/suppl/GSE142241_RAW.tar;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121721/suppl/GSE121721_RAW.tar;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121720/suppl/GSE121720%5FRNAseq%5Fexpression%5Fmatrix%5FTPMs.txt.gz;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/matrix/GSE186458_series_matrix.txt.gz;

gzip -d GSE222147_lcm_wbgs_bsseq_smoothed_coverage_filtered_geo.rda.gz;

mkdir GSE121721_RAW; cd GSE121721_RAW;
mv ../GSE121721_RAW.tar ./;
tar -xvf GSE121721_RAW.tar;
for i in `ls *methylation_values.bigWig`; do bigWigToBedGraph $i $i.bedGraph; done
cd ..;

mkdir GSE142241_RAW; cd GSE142241_RAW;
mv ../GSE142241_RAW.tar ./;
tar -xvf GSE142241_RAW.tar;
for i in `ls *.bw`; do bigWigToBedGraph $i $i.bedGraph; done
cd ..;

for i in `zcat GSE186458_series_matrix.txt.gz |grep "Sample_supplementary_file_1"`; do echo $i >> GSE186458.list; done
mkdir GSE186458; cd GSE186458;
for i in `grep beta ../GSE186458.list|sed 's/ftp:/https:/g' |cut -f2 -d '"'`; do wget $i; done
for i in *.beta; do wgbstools view $i|awk '$5>=10 && $5 <=150' OFS="\t" $F | perl -ane '$F[5]=$F[3]/$F[4]; print "$F[0]\t$F[1]\t$F[2]\t$F[5]\n"' | grep -v chrX | grep -v chrY | grep -v chrM > $i.bedGraph; done
cd ..;

wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz
wget https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA_GTEX_category.txt

cd ../..; python ./preprocess.py
rm data/geo/*/*.bedGraph
rm data/geo/*/*.beta
rm data/geo/*/*.bw
rm data/geo/*/*.bigWig
rm data/**/*.m
rm data/**/*.cov
rm data/**/*.pos
rm data/**/*.tab
rm data/geo/*/*.tar
rm data/geo/*.rda
