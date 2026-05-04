mkdir data;cd data;
mkdir geo;cd geo;

wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121720/suppl/GSE121720_RAW.tar;

mkdir GSE121720_RAW; cd GSE121720_RAW;
mv ../GSE121720_RAW.tar ./;
tar -xvf GSE121720_RAW.tar;
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz;
gzip -d gencode.v19.annotation.gtf.gz;
gtf2bed < gencode.v19.annotation.gtf | cut -f1-8 |grep gene > gencode.v19.annotation.onlygenes.bed;

multiBigwigSummary BED-file -p 16 -v -b *.bigWig -o ../results.genes.npz --BED gencode.v19.annotation.onlygenes.bed --outRawCounts ../raw.genes.npz
