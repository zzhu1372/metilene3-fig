# %%
import os
import pandas as pd
from numpy import random
import seaborn as sns
import matplotlib.pyplot as plt

# %%
wd = os.getcwd()
# os.system('tar -xzvf bg.tar.gz')
wd

# %%
def scale(x,a):
    return (x-0.5)*a+0.5

def min0max1(x):
    return min(1,max(0,x))
    
def dmr_simulate(bgab, abab, rscript, a, bg):
    ### simulate foreground data
    os.system('mkdir -p ./data/'+bgab.replace(' ','')+'-fg-a'+str(a))
    os.system('Rscript '+wd+'/'+rscript+' '+\
                        abab+' '+bg+' '+\
                          wd+'/data/'+bgab.replace(' ','')+'-fg-a'+str(a)+' '+str(a)+'')


    
    ### read DMRs and foreground data
    dmrs = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a10/DMRs_'+abab.replace(' ','_')+'_a10.bed',header=None)
    beta = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a'+str(a)+'/beta_'+abab.replace(' ','_')+'_a'+str(a)+'.tsv')

    
    
    ### simulate LGGGBM purity
    ### https://doi.org/10.1038/ncomms9971
    tcgaPurity = pd.read_excel('https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx',skiprows=3)
    tcgaPurity = tcgaPurity.loc[tcgaPurity['Cancer type'].isin(['LGG','GBM'])] 
    purityMean = tcgaPurity['CPE'].dropna().mean()
    purityStd = tcgaPurity['CPE'].dropna().std()
    print(purityMean, purityStd)
    random.seed(0)
    purity = random.normal(loc=purityMean, scale=purityStd, size=(1, 50))[0]
    for i in range(len(purity)):
        purity[i] = max(0,purity[i])
        purity[i] = min(1,purity[i])
    print(min(purity), max(purity))
    

    
    ### simulate methylation of C1 C2 regions
    nogrpdmrs = dmrs.loc[dmrs[6]=='nogrp']
    random.seed(1)
    c1factor = random.uniform(0,1,50)
    c2factor = random.normal(0.5,0.1,50)
    for i in beta.columns[2:]:
        j = int(i.split('s')[-1])
        random.seed(j)
        nogrpdmrs[i] = \
            sum([\
                sum([[min0max1(x) for x in random.normal(loc=scale(c1factor[j],0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(1-c1factor[j],0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(min0max1(c2factor[j]),0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(1-min0max1(c2factor[j]),0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[])\
                ],[])
    nogrpdmrs.index = nogrpdmrs[4]
    nogrpdmrs.index.name = None
    
    ### replace CpG met
    cpgs = beta[['chr','pos']]
    cpgs['start'] = cpgs['pos']-1
    
    from pybedtools import BedTool
    def find_overlapping_regions_df(bed1_df, bed2_df, bed1_cols, bed2_cols):
        bed_1 = BedTool.from_dataframe(bed1_df[bed1_cols].sort_values(bed1_cols))
        bed_2 = BedTool.from_dataframe(bed2_df[bed2_cols].sort_values(bed2_cols))
        return BedTool.to_dataframe(bed_1.intersect(bed_2, wa=True, wb=True))
    
    nogrpcpgs = find_overlapping_regions_df(nogrpdmrs, cpgs, \
                                            [0,1,2,4], ['chr','start','pos'])
    print((nogrpdmrs[5]-nogrpdmrs[4]+1).sum())
    
    for i in nogrpdmrs.columns[nogrpdmrs.columns.astype(str).str.contains('_s')]:
        nogrpcpgs[i] = nogrpcpgs['name'].map(nogrpdmrs[i])
    
    print(beta.loc[beta['pos'].isin(nogrpcpgs['thickStart'])].shape)
    nogrpcpgs.index = nogrpcpgs['thickStart']
    
    finalfactor = beta[[]]
    finalfactor.index = beta['pos']

    ### merge purity and C1 C2 regions
    for i in nogrpdmrs.columns[nogrpdmrs.columns.astype(str).str.contains('_s')]:
        j = int(i.split('s')[-1])
        finalfactor[i] = purity[j]
    
    finalfactor = finalfactor.loc[~finalfactor.index.isin(set(nogrpcpgs.index))]
    finalfactor = pd.concat([finalfactor,nogrpcpgs[finalfactor.columns]]).sort_index()



    ### simulate two batches
    met = beta[['chr','pos']]
    bcpgs0 = set(met.sample(500, random_state=0).sort_index()['pos'])
    bcpgs1 = set(met.sample(500, random_state=1).sort_index()['pos'])
    
    

    ### generate final methylation rates
    ctrs = sorted(pd.Series(os.listdir(bg)).loc[\
        pd.Series(os.listdir(bg)).str.contains('sample')])
    
    purityDict = {}
    c1factorDict = {}
    c2factorDict = {}

    os.mkdir('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a))
    for i in beta.columns[2:]:
        j = int(i.split('s')[-1])
        tmp = pd.read_table(bg+ctrs[j],header=None, index_col=1)
        bed = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a'+str(a)+'/'+i+'.beta_'+abab.replace(' ','_')+'.bed', header=None)
        purityDict[i] = purity[j]
        c1factorDict[i] = c1factor[j]
        c2factorDict[i] = c2factor[j]
    
        if j%2==0:
            bed[3] = (bed[3]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[4])*(1-bed[1].map(finalfactor[i]))).astype(int) * (1-bed[1].isin(bcpgs0))
            # print(i,j,'b1')
        else:
            bed[3] = (bed[3]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[4])*(1-bed[1].map(finalfactor[i]))).astype(int) * (1-bed[1].isin(bcpgs1))
            # print(i,j)
            
        bed[4] = (bed[4]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[5])*(1-bed[1].map(finalfactor[i]))).astype(int)
        bed[5] = bed[3]/bed[4]
        bed.index = bed[1]
        met[i] = met['pos'].map(bed[5])
        wgbstools = bed[[0,1,2,3,4]]
        methylscore = bed[[0,1,2,5,3,4]]
        methylscore[4] = methylscore[4] - methylscore[3]
        methylscore[1] = methylscore[1] - 1
        methylscore[2] = methylscore[2] - 1
        wgbstools.to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'a'+str(a)+'.'+i+'.wgbstools.bed', sep='\t', float_format='%.3f', header=False, index=False)
        methylscore.to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'a'+str(a)+'.'+i+'.methylscore.bed', sep='\t', float_format='%.3f', header=False, index=False)
    
    pd.Series(purityDict).to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.purity.csv')
    pd.Series(c1factorDict).to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.c1.csv')
    pd.Series(c2factorDict).to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.c2.csv')
    
    met.to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'a'+str(a)+'.metilene3.tsv', sep='\t', float_format='%.3f', na_rep='.', index=False)

    ### subsample for consensus analysis
    if bgab.replace(' ','')=='403':
        os.mkdir('./data/'+bgab.replace(' ','')+'-consensus')
        for i in range(50):
            asize = 1+i%10
            selected = sorted(pd.Series(met.columns[met.columns.str.contains('0_')]).sample(asize, random_state=i))+\
                        sorted(pd.Series(met.columns[~met.columns.str.contains('0_|chr|pos')]).sample(40-asize, random_state=i))
            met[['chr','pos']+selected].to_csv('./data/'+bgab.replace(' ','')+'-consensus/'+bgab.replace(' ','')+'a'+str(a)+'.metilene3.n40rep'+str(i)+'.tsv',\
                   sep='\t', float_format='%.3f', na_rep='.', index=False)
    
    met['End'] = met['pos']+1
    met = met[['chr','pos','End']+sorted(met.columns[met.columns.str.contains('_')])]
    newids = []
    j = 0
    g = '0'
    for i in sorted(met.columns[met.columns.str.contains('_')]):
        if i.split('_')[0]==g:
            j+=1
        else:
            g = i.split('_')[0]
            j=1
        newids.append('G'+i.split('_')[0]+'_'+str(j))
    
    met.columns = ['Chrome','Start','End']+newids
    met.to_csv('./data/'+bgab.replace(' ','')+'-LGGGBM-a'+str(a)+'/'+bgab.replace(' ','')+'a'+str(a)+'.smart2.tsv', sep='\t', float_format='%.3f', na_rep='.', index=False)

# %%
bgab = '40 3'
abab = bgab+' 22 22'
rscript = './simulate-grpmet.R'
a = 10

dmr_simulate(bgab, abab, rscript, a, './bg/'+''.join(bgab.split(' '))+'/')

# %%
bgab = '15 5'
abab = bgab+' 10 10'
rscript = './simulate-grpmet.R'
a = 10

dmr_simulate(bgab, abab, rscript, a, './bg/'+''.join(bgab.split(' '))+'/')

# %%
def dmr_lowpurity_simulate(bgab, abab, rscript, a, bg, purityMean):
    ### read DMRs and foreground data
    dmrs = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a10/DMRs_'+abab.replace(' ','_')+'_a10.bed',header=None)
    beta = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a'+str(a)+'/beta_'+abab.replace(' ','_')+'_a'+str(a)+'.tsv')

    
    
    ### simulate LGGGBM purity
    purityStd = 0.1
    print(purityMean, purityStd)
    random.seed(0)
    purity = random.normal(loc=purityMean, scale=purityStd, size=(1, 50))[0]
    for i in range(len(purity)):
        purity[i] = max(0,purity[i])
        purity[i] = min(1,purity[i])
    print(min(purity), max(purity))
    

    
    ### simulate methylation of C1 C2 regions
    nogrpdmrs = dmrs.loc[dmrs[6]=='nogrp']
    random.seed(1)
    c1factor = random.uniform(0,1,50)
    c2factor = random.normal(0.5,0.1,50)
    for i in beta.columns[2:]:
        j = int(i.split('s')[-1])
        random.seed(j)
        nogrpdmrs[i] = \
            sum([\
                sum([[min0max1(x) for x in random.normal(loc=scale(c1factor[j],0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(1-c1factor[j],0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(min0max1(c2factor[j]),0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[]),\
                sum([[min0max1(x) for x in random.normal(loc=scale(1-min0max1(c2factor[j]),0.1*(pho+1)), scale=0.3, size=(1, int(len(nogrpdmrs.index)/40)))[0]] for pho in range(10)],[])\
                ],[])
    nogrpdmrs.index = nogrpdmrs[4]
    nogrpdmrs.index.name = None
    
    ### replace CpG met
    cpgs = beta[['chr','pos']]
    cpgs['start'] = cpgs['pos']-1
    
    from pybedtools import BedTool
    def find_overlapping_regions_df(bed1_df, bed2_df, bed1_cols, bed2_cols):
        bed_1 = BedTool.from_dataframe(bed1_df[bed1_cols].sort_values(bed1_cols))
        bed_2 = BedTool.from_dataframe(bed2_df[bed2_cols].sort_values(bed2_cols))
        return BedTool.to_dataframe(bed_1.intersect(bed_2, wa=True, wb=True))
    
    nogrpcpgs = find_overlapping_regions_df(nogrpdmrs, cpgs, \
                                            [0,1,2,4], ['chr','start','pos'])
    print((nogrpdmrs[5]-nogrpdmrs[4]+1).sum())
    
    for i in nogrpdmrs.columns[nogrpdmrs.columns.astype(str).str.contains('_s')]:
        nogrpcpgs[i] = nogrpcpgs['name'].map(nogrpdmrs[i])
    
    print(beta.loc[beta['pos'].isin(nogrpcpgs['thickStart'])].shape)
    nogrpcpgs.index = nogrpcpgs['thickStart']
    
    finalfactor = beta[[]]
    finalfactor.index = beta['pos']

    ### merge purity and C1 C2 regions
    for i in nogrpdmrs.columns[nogrpdmrs.columns.astype(str).str.contains('_s')]:
        j = int(i.split('s')[-1])
        finalfactor[i] = purity[j]
    
    finalfactor = finalfactor.loc[~finalfactor.index.isin(set(nogrpcpgs.index))]
    finalfactor = pd.concat([finalfactor,nogrpcpgs[finalfactor.columns]]).sort_index()



    ### simulate two batches
    met = beta[['chr','pos']]
    bcpgs0 = set(met.sample(500, random_state=0).sort_index()['pos'])
    bcpgs1 = set(met.sample(500, random_state=1).sort_index()['pos'])
    
    

    ### generate final methylation rates
    ctrs = sorted(pd.Series(os.listdir(bg)).loc[\
        pd.Series(os.listdir(bg)).str.contains('sample')])
    
    purityDict = {}
    c1factorDict = {}
    c2factorDict = {}

    os.mkdir('./data/'+bgab.replace(' ','')+'-puritymean'+str(int(purityMean*100))+'-a'+str(a))
    for i in beta.columns[2:]:
        j = int(i.split('s')[-1])
        tmp = pd.read_table(bg+ctrs[j],header=None, index_col=1)
        bed = pd.read_table('./data/'+bgab.replace(' ','')+'-fg-a'+str(a)+'/'+i+'.beta_'+abab.replace(' ','_')+'.bed', header=None)
        purityDict[i] = purity[j]
        c1factorDict[i] = c1factor[j]
        c2factorDict[i] = c2factor[j]
    
        if j%2==0:
            bed[3] = (bed[3]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[4])*(1-bed[1].map(finalfactor[i]))).astype(int) * (1-bed[1].isin(bcpgs0))
            # print(i,j,'b1')
        else:
            bed[3] = (bed[3]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[4])*(1-bed[1].map(finalfactor[i]))).astype(int) * (1-bed[1].isin(bcpgs1))
            # print(i,j)
            
        bed[4] = (bed[4]*bed[1].map(finalfactor[i]) + bed[1].map(tmp[5])*(1-bed[1].map(finalfactor[i]))).astype(int)
        bed[5] = bed[3]/bed[4]
        bed.index = bed[1]
        met[i] = met['pos'].map(bed[5])
        
    pd.Series(purityDict).to_csv('./data/'+bgab.replace(' ','')+'-puritymean'+str(int(purityMean*100))+'-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.purity.csv')
    pd.Series(c1factorDict).to_csv('./data/'+bgab.replace(' ','')+'-puritymean'+str(int(purityMean*100))+'-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.c1.csv')
    pd.Series(c2factorDict).to_csv('./data/'+bgab.replace(' ','')+'-puritymean'+str(int(purityMean*100))+'-a'+str(a)+'/'+bgab.replace(' ','')+'.a'+str(a)+'.c2.csv')
    
    met.to_csv('./data/'+bgab.replace(' ','')+'-puritymean'+str(int(purityMean*100))+'-a'+str(a)+'/'+bgab.replace(' ','')+'a'+str(a)+'.metilene3.tsv', sep='\t', float_format='%.3f', na_rep='.', index=False)

# %%
bgab = '40 3'
abab = bgab+' 22 22'
rscript = './simulate-grpmet.R'
a = 10

for purityMean in [0.5,0.6,0.7,0.8]:
    print(purityMean)
    dmr_lowpurity_simulate(bgab, abab, rscript, a, './bg/'+''.join(bgab.split(' '))+'/', purityMean)

# %%
bgab = '15 5'
abab = bgab+' 10 10'
rscript = './simulate-grpmet.R'
a = 10

for purityMean in [0.5,0.6,0.7,0.8]:
    print(purityMean)
    dmr_lowpurity_simulate(bgab, abab, rscript, a, './bg/'+''.join(bgab.split(' '))+'/', purityMean)
