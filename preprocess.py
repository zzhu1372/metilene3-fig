import pandas as pd
import random
import os


# Loyfer
# https://doi.org/10.1038/s41586-022-05580-6
df = pd.read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM4_ESM.xlsx', skiprows=2)

bedfiles = {}
for i in df['Sample name']:
    x = i
    for j in os.listdir('./data/geo/GSE186458/'):
        if j.split('.')[-1]!='bedGraph':
            continue
        if j.find(i.split('-')[-1])!=-1:
            x = j
            break
    if x==i:
        x = None
    else:
        bedfiles[i] = x

df['bedfile'] = df['Sample name'].map(bedfiles)
df = df.sort_values(['Group','Sample name'])
df.index = range(len(df.index))

random.seed(1)
new_names = list(df.index)
random.shuffle(new_names)
df['ID'] = ['S'+str(i) for i in new_names]
df['ID'] = df['ID'] +'='+ df['Sample name']
df = df.sort_values('ID')

os.system('cd ./data/geo/GSE186458/;\
bedtools unionbedg -i '+(df['bedfile']+' ').sum()+\
'-header -names '+(df['ID']+' ').sum()+\
'-filler .| '+os.getcwd()+'/round3.awk |cut -f 1,3- > '+os.getcwd()+'/data/GSE186458_celltypes.input.tsv')

df = df.loc[df['ID'].str.contains('Blood')]

os.system('cd ./data/geo/GSE186458/;\
bedtools unionbedg -i '+(df['bedfile']+' ').sum()+\
'-header -names '+(df['ID']+' ').sum()+\
'-filler .| '+os.getcwd()+'/round3.awk |cut -f 1,3- > '+os.getcwd()+'/data/GSE186458_blood.input.tsv')


# Wu
# https://doi.org/10.1038/s41467-020-20225-w
df = pd.read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-20225-w/MediaObjects/41467_2020_20225_MOESM3_ESM.xlsx')
df['IDH1 mutaton status'] = df['IDH1 mutaton status'].fillna('normal')

bedfiles = {}

normal = {
    'NBr1':'frontal',
    'NBr2':'occipital',
    'NBr3':'parietal',
    'NBr4':'temporal'
}
for i in df['Sample ID']:
    x = i
    for j in os.listdir('./data/geo/GSE121721_RAW/'):
        if i.find('NBr')==-1:
            if j.find(i+'_methylation')!=-1 and j.find('.bedGraph')!=-1:
                x = j
                break
        if i.find('NBr')!=-1:
            if j.find(normal[i])!=-1 and j.find('.bedGraph')!=-1:
                x = j
                break
    if x==i:
        print('err')
    else:
        bedfiles[i] = x

df['bedfile'] = df['Sample ID'].map(bedfiles)
df = df.sort_values(['IDH1 mutaton status','Sample ID'])
df.index = range(len(df.index))

random.seed(1)
new_names = list(df.index)
random.shuffle(new_names)
df['ID'] = ['S'+str(i) for i in new_names]
df['ID'] = df['ID'] +'='+ df['IDH1 mutaton status']+'_'+df['Sample ID']
df = df.sort_values('ID')

os.system('cd ./data/geo/GSE121721_RAW/;\
bedtools unionbedg -i '+(df['bedfile']+' ').sum()+\
'-header -names '+(df['ID']+' ').sum()+\
'-filler .| '+os.getcwd()+'/round3.awk |cut -f 1,3- > '+os.getcwd()+'/data/GSE121721_glioma.input.tsv')


# Li
flist = pd.Series(os.listdir('./data/geo/GSE142241_RAW/'))
flist = flist.loc[flist.str.contains('.bw.bedGraph')&flist.str.contains('CSF')\
                    &(~flist.str.contains('CMS'))]
rename = {'chrom':'chrom','end':'end'}
random.seed(0)
random_num = list(range(len(flist)))
random.shuffle(random_num)
for ii,i in enumerate(flist):
    rename[i] = 'S'+str(random_num[ii])+'='+'_'.join(i.split('.')[0].split('_')[1:])

flist = pd.DataFrame([rename[i] for i in flist], flist).sort_values(0)
os.system('cd ./data/geo/GSE142241_RAW/;\
bedtools unionbedg -i '+' '.join(list(flist.index))+\
' -header -names '+' '.join(list(flist[0]))+\
' -filler .| '+os.getcwd()+'/round3.awk |cut -f 1,3- > '+os.getcwd()+'/data/GSE142241_cfMB.input.tsv')


# Lo
os.system("Rscript -e \"library('bsseq');\
load('./data/geo/GSE222147_lcm_wbgs_bsseq_smoothed_coverage_filtered_geo.rda');\
write.table(getMeth(bsseq_smoothed_coverage_filtered, type = 'raw', what = 'perBase'),'./data/geo/GSE222147.tab', \
row.names = FALSE, quote = FALSE, sep = '\\t');\
write.table(granges(bsseq_smoothed_coverage_filtered),'./data/geo/GSE222147.pos', \
row.names = FALSE, quote = FALSE, sep = '\\t');\
write.table(getCoverage(bsseq_smoothed_coverage_filtered),'./data/geo/GSE222147.cov', \
row.names = FALSE, quote = FALSE, sep = '\\t');\
write.table(getCoverage(bsseq_smoothed_coverage_filtered, type = 'M'),'./data/geo/GSE222147.m', \
row.names = FALSE, quote = FALSE, sep = '\\t');\"")

met = pd.concat([pd.read_table('./data/geo/GSE222147.pos')[['seqnames','end']],\
       pd.read_table('./data/geo/GSE222147.tab')], axis=1,)

met['seqnames'] = 'chr'+met['seqnames'].astype(str)
met = met.sort_values(['seqnames','end'])

new_names = list(met.columns[2:])
random.seed(1)
random.shuffle(new_names)
met = met[list(met.columns[:2])+new_names]
for i in range(len(new_names)):
    new_names[i] = 'S'+str(i)+'='+new_names[i]
met.columns = ['chr','pos']+new_names
met = met[list(met.columns[:2])+sorted(new_names)]

met.to_csv('./data/GSE222147_pdac.input.tsv', sep='\t', float_format='%.3f', na_rep='.', index=False)
