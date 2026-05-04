# metilene3-fig
Codes to generate figures.
```
git clone https://github.com/zzhu1372/metilene3-fig.git
cd metilene3-fig
mkdir SourceData
figwd=$PWD
```


## Simulation
1. Create folders
```
cd $figwd/simulation
mkdir figures
mkdir conda-envs
cd conda-envs
envswd=$PWD
```

2. Simulate two backgrounds
```
cd $envswd
conda create --prefix ./bgsimulation-env -y -c conda-forge -c bioconda python==3.10.0 r-base==4.4.3
conda activate ./bgsimulation-env
cd ..
# download the chr10 annotation
wget https://genome.cshlp.org/content/suppl/2015/12/07/gr.196394.115.DC1/metilene_DMR_simulation_scripts.tar.gz
tar -xzvf metilene_DMR_simulation_scripts.tar.gz
mv metilene_DMR_simulation_scripts/chromatin_annotation_chr10.txt ./
mkdir -p bg/403
for i in {0..49}; do echo $i;Rscript $PWD/simulate_bg.R 40 3 $PWD/chromatin_annotation_chr10.txt $i > $PWD/bg/403/sample_$i.txt; done
mkdir -p bg/155
for i in {0..49}; do echo $i;Rscript $PWD/simulate_bg.R 15 5 $PWD/chromatin_annotation_chr10.txt $i > $PWD/bg/155/sample_$i.txt; done
tar -czvf bg.tar.gz ./bg
```

3. Simulate methylation rates
```
cd $envswd
conda create --prefix ./simulation-env -y -c conda-forge -c bioconda python==3.10.0 pandas==1.5.1 openpyxl==3.1.5 scikit-learn==1.3.1 seaborn==0.12.1 biopython==1.85 pybedtools==0.12.0 r-base==4.4.3 r-hash==3.0.1 r-plyr==1.8.9 jupyterlab==4.4.5
conda activate ./simulation-env
cd ..
python 0-simulate-DMR.py
tar -czvf simulated-data.tar.gz ./data
```

4. DMR identification
```
cd $envswd
conda create --prefix ./smart2-env -y -c conda-forge -c bioconda python==2.7.15 pip==20.1.1 numpy==1.16.0 scipy==1.2.0 statsmodels==0.8.0
conda activate ./smart2-env
pip install SMART-BS-Seq==2.2.8

cd $envswd
conda create --prefix ./wgbstools-env -y -c conda-forge -c bioconda python==3.10.0 pandas==1.5.1 numpy==1.26.4 scipy==1.15.2 samtools==1.23.1 tabix==1.11
conda activate ./wgbstools-env
cd $(which python| xargs -0 dirname)
wget https://github.com/nloyfer/wgbs_tools/archive/refs/tags/0.2.2.zip
unzip 0.2.2.zip
cd wgbs_tools-0.2.2
python setup.py
ln -s $(realpath ./wgbstools) ../wgbstools
wgbstools init_genome hg19 -@ 1

cd $envswd
conda create --prefix ./methylscore-env -y -c conda-forge -c bioconda python==3.10.0 nextflow==25.10.4
conda activate ./methylscore-env
cd $(which python| xargs -0 dirname)
wget https://github.com/Computomics/MethylScore/archive/refs/tags/0.2.zip
unzip 0.2.zip
cd MethylScore-0.2
singularity pull methylscore.sif docker://quay.io/beckerlab/methylscore@sha256:2de46d1fdf9dd48770121d25ff3224fea64c7fddc25aa9edca3bcad675156754
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/chr10.fa.gz
gzip -d chr10.fa.gz

cd $envswd
conda create --prefix ./metilene3-env -y -c conda-forge -c bioconda python==3.10.0 pandas==1.5.1 scipy==1.15.2
conda activate ./metilene3-env
cd $figwd/simulation
git clone https://github.com/zzhu1372/metilene3.git
cd metilene3
make
cd ..
```
Adjust the ``conda_init`` variable and the HPC job submission command in ``1-dmrs-run.py`` with your own configuration.
```
python 1-dmrs-run.py
```

5. Analysis
```
cd $envswd
conda activate ./simulation-env
cd ..
```
Run:
```
0-plot-DMR.ipynb
1-dmrs.ipynb
2-dmtree.ipynb
```


## Biological datasets
1. Download conda environment and public data
```
cd $figwd/biological
mkdir figures
mkdir conda-envs
cd conda-envs
envswd=$PWD

cd $envswd
conda create --prefix ./deeptools-env -y -c conda-forge -c bioconda deeptools==3.5.6
cd ./deeptools-env
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
cd $envswd
conda create --prefix ./biological-env -y -c conda-forge -c bioconda perl==5.32.1 python==3.10.0 pandas==1.5.1 openpyxl==3.1.5 scikit-learn==1.3.1 seaborn==0.12.1 biopython==1.85 pybedtools==0.12.0 gseapy==1.1.9 r-base==4.4.3 r-hash==3.0.1 r-plyr==1.8.9 r-rcolorbrewer==1.1_3 bioconductor-enrichedheatmap==1.36.0 r-circlize==0.4.16 bioconductor-bsseq==1.42.0 bioconductor-deseq2==1.46.0 bioconductor-ChIPseeker==1.42.0 bioconductor-org.Hs.eg.db==3.20.0 bioconductor-txdb.hsapiens.ucsc.hg19.knowngene==3.2.2 homer==5.1 samtools==1.23.1 tabix==1.11 jupyterlab==4.4.5
conda activate ./biological-env
cd ./biological-env
wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.41.tar.bz2
cd bin
wget https://hgdownload.gi.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod +x bigWigToBedGraph 
wget https://github.com/nloyfer/wgbs_tools/archive/refs/tags/0.2.2.zip
unzip 0.2.2.zip
cd wgbs_tools-0.2.2
python setup.py
ln -s $(realpath ./wgbstools) ../wgbstools
cd ../../../..

wgbstools init_genome hg19 -@ 1 -f

perl $(dirname $(which findMotifs.pl))/../share/homer/configureHomer.pl -install hg19

git clone https://github.com/zzhu1372/metilene3.git
cd metilene3
sed 's/-O3/-O0/g' Makefile > Makefile.1
mv Makefile.1 Makefile
make
cd ..

cd $figwd/biological
chmod +x round3.awk
chmod +x get_data_step1.sh
chmod +x get_data_step2.sh
conda activate ./conda-envs/deeptools-env
./get_data_step1.sh
conda activate ./conda-envs/biological-env
./get_data_step2.sh
tar -czvf biological-data.tar.gz ./data
```

2. Analysis

Run:
```
analysis-blood.ipynb
analysis-gbm.ipynb
analysis-pdac.ipynb
analysis-cfMB.ipynb
analysis-celltypes.ipynb
```
