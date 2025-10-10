# metilene3-fig
Codes to generate figures.

## Dependencies
1. metilene<sup>3</sup> (should be installed under the same folder)
1. wgbstools (version 0.2.0)
1. SMART (version 2.2.8)
1. scipy (version 1.15.2)
1. datatable (version 1.1.0)
1. bedtools (version 2.31.1)
1. pybedtools (version 0.12.0)
1. HOMER (version 5.1)
1. DESeq2 (version 1.46.0)
1. bsseq (version 1.42.0)

## Download public data
```
git clone https://github.com/zzhu1372/metilene3-fig.git
cd metilene3-fig
chmod +x round3.awk
bash get_data.sh
```

## Run the analysis

``analysis-simulation.ipynb``: related to figure 2 (run ``simulation/simulations.ipynb`` and then ``simulation/benchmarking.ipynb`` before to get the results on the simulated dataset.)

``analysis-blood.ipynb``: related to figure 3

``analysis-gbm.ipynb``: related to figure 4

``analysis-pdac.ipynb``: related to figure 5 and 6

