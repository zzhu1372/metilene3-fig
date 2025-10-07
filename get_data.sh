mkdir data;cd data;
mkdir geo;cd geo;

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE222147/suppl/GSE222147_lcm_wbgs_bsseq_smoothed_coverage_filtered_geo.rda.gz;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE142nnn/GSE142241/suppl/GSE142241_RAW.tar;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121721/suppl/GSE121721_RAW.tar;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121720/suppl/GSE121720_RAW.tar;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121720/suppl/GSE121720%5FRNAseq%5Fexpression%5Fmatrix%5FTPMs.txt.gz;
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186458/matrix/GSE186458_series_matrix.txt.gz;

gzip -d GSE222147_lcm_wbgs_bsseq_smoothed_coverage_filtered_geo.rda.gz;

mkdir GSE121721_RAW; cd GSE121721_RAW;
mv ../GSE121721_RAW.tar ./;
tar -xvf GSE121721_RAW.tar;
for i in `ls *methylation_values.bigWig`; do bigWigToBedGraph $i $i.bedGraph; done
cd ..;

mkdir GSE121720_RAW; cd GSE121720_RAW;
mv ../GSE121720_RAW.tar ./;
tar -xvf GSE121720_RAW.tar;
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz;
gzip -d gencode.v19.annotation.gtf.gz;
gtf2bed < gencode.v19.annotation.gtf | cut -f1-8 |grep gene > gencode.v19.annotation.onlygenes.bed;
multiBigwigSummary BED-file -p 16 -v -b *.bigWig -o ../results.genes.npz --BED gencode.v19.annotation.onlygenes.bed --outRawCounts ../raw.genes.npz
cd ..;

mkdir GSE142241_RAW; cd GSE142241_RAW;
mv ../GSE142241_RAW.tar ./;
tar -xvf GSE142241_RAW.tar;
for i in `ls *.bw`; do bigWigToBedGraph $i $i.bedGraph; done
cd ..;

for i in `zcat GSE186458_series_matrix.txt.gz |grep "Sample_supplementary_file_1"`; do echo $i >> GSE186458.list; done
mkdir GSE186458; cd GSE186458;
for i in `grep Blood ../GSE186458.list |cut -f2 -d '"'`; do wget $i; done
wgbstools init_genome hg19 -f
for i in *.beta; do wgbstools view $i|awk '$5>=10 && $5 <=150' OFS="\t" $F | perl -ane '$F[5]=$F[3]/$F[4]; print "$F[0]\t$F[1]\t$F[2]\t$F[5]\n"' | grep -v chrX | grep -v chrY | grep -v chrM > $i.bedGraph; done
cd ..;

cd ../..; python ./preprocess.py