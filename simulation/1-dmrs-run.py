# %%
import os
import pandas as pd

conda_init = '''
cd
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/zzhu/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/zzhu/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/zzhu/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/zzhu/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


# >>> mamba initialize >>>
# !! Contents within this block are managed by 'mamba shell init' !!
export MAMBA_EXE='/home/zzhu/miniforge3/bin/mamba';
export MAMBA_ROOT_PREFIX='/home/zzhu/miniforge3';
__mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"  # Fallback on help from mamba activate
fi
unset __mamba_setup
# <<< mamba initialize <<<

conda activate ~/.conda/envs/base
'''

# %%
wd = os.getcwd()
wgbstools = wd+'/conda-envs/wgbstools-env'
metilene = wd+'/conda-envs/metilene3-env'
ms = wd+'/conda-envs/methylscore-env'
smart = wd+'/conda-envs/smart2-env'

# %%
def run_all(ab):
    ab = ''.join(ab.split('_')[:2])
    indir = wd+'/data/'+ab+'-LGGGBM-a10/'

    s = pd.DataFrame(os.listdir(indir))
    s = s.loc[s[0].str.contains('wgbstools.bed')&s[0].str.contains(ab+'a10')]
    s[0] = s[0].apply(lambda x:x.replace('.bed',''))
    s[1] = s[0].apply(lambda x:'G'+x.split('_')[0].split('.')[-1])
    s.columns = ['name','group']
    s.to_csv(indir+'/'+ab+'.groups.wgbstools.csv',index=False)
    s.columns = ['ID','Group']
    s['ID'] = s['ID'].apply(lambda x:x.split('.')[1])
    s.to_csv(indir+'/'+ab+'.groups.metilene.tsv',index=False,sep='\t')
    outdir = wd+'/data/results/sup-'+ab

    os.system('mkdir -p '+outdir)
    
    mxqsh = outdir+"/"+ab+".mxq.sh"
    
    with open(mxqsh,'w') as fn:
        fn.write(conda_init)
        fn.write("cd $MXQ_JOB_TMPDIR/;mkdir input;cp -r "+indir+"/* input/;\n")
        fn.write('cp '+ms+'/bin/MethylScore-0.2/chr10.fa ./;\n')
        fn.write('conda activate '+wgbstools+';\n')
        fn.write('wgbstools bed2beta input/'+ab+'a10.*.wgbstools.bed;\n')
        for i in [0]:
            for c in ['10','1']:
                fn.write("mkdir rep"+str(i)+'.c'+c+";cd rep"+str(i)+'.c'+c+";\n")
                fn.write('''
python -c "import os
import pandas as pd
samples = []
beds = []
for i in os.listdir('../input/'):
    if i.split('.')[-1]!='bed': continue
    if i.split('.')[-2]!='methylscore': continue
    samples.append(i.split('.')[1])
    beds.append(os.getcwd()+'/../input/'+i)
df = pd.DataFrame(beds,samples)
df.to_csv('./sample.txt',sep='\t',header=False)" \n
''')
                
                fn.write("conda activate "+ms+";\n")
                fn.write('''nextflow run '''+ms+'''/bin/MethylScore-0.2/main.nf \
--BEDGRAPH --SAMPLE_SHEET='''+'$PWD/sample.txt'+'''  \
--GENOME=$PWD/../chr10.fa  --DMR_CONTEXTS=CG \
--PROJECT_FOLDER=$PWD/ms-prep/ \
-w '''+'$PWD/ms-prep-work \n') 
                fn.write('mv ./ms-prep-work/tmp/*/*/genome_matrix.tsv ./;\n')
                fn.write('''/usr/bin/time -v -o timeNmem.ms nextflow run '''+ms+'''/bin/MethylScore-0.2/main.nf \
-with-singularity '''+ms+'''/bin/MethylScore-0.2/methylscore.sif  \
--MATRIX='''+'$PWD/genome_matrix.tsv \
--GENOME=$PWD/../chr10.fa  --DMR_CONTEXTS=CG \
--PROJECT_FOLDER='''+'$PWD/ms-results/ \
-w '''+'$PWD/ms-results-work \n\n\n') 
                fn.write('rm -rf ./*work \n')
                fn.write("cd ..;cp -rf rep"+str(i)+'.c'+c+" "+outdir+"/ \n")
                fn.write("cd rep"+str(i)+'.c'+c+" \n\n\n")

                
                fn.write("conda activate "+wgbstools+";\n")
                fn.write('/usr/bin/time -v -o timeNmem.wgbstools_seg wgbstools segment --betas ../'+ab+'a10.*.wgbstools.beta -@ '+c+' -r chr10 -o ./blocks.small.bed\n')
                fn.write('/usr/bin/time -v -o timeNmem.wgbstools_dmr wgbstools find_markers -@ '+c+' --blocks_path '+\
                          './blocks.small.bed --groups_file '+\
                          '../input/'+ab+'.groups.wgbstools.csv --betas '+\
                          '../'+ab+'a10.*.wgbstools.beta -o ./wgbstools \n')
                fn.write("cd ..;cp -rf rep"+str(i)+'.c'+c+" "+outdir+"/ \n")
                fn.write("cd rep"+str(i)+'.c'+c+" \n\n\n")

                
                fn.write("conda activate "+smart+";\n")
                fn.write("/usr/bin/time -v -o timeNmem.smart SMART ../input/"+ab+"a10.smart2.tsv \
                -t DeNovoDMR \
                -o ./smart \n")
                
                
                fn.write("conda activate "+metilene+";\n")
                fn.write('/usr/bin/time -v -o timeNmem.metilene_sup python '+wd+'/metilene3/metilene3.py \
                        -i ../input/'+ab+'a10.metilene3.tsv \
                        -g ../input/'+ab+'.groups.metilene.tsv \
                        -o ./metilene_sup \
                        -t '+c+'\n')
                fn.write("cd ..;cp -rf rep"+str(i)+'.c'+c+" "+outdir+"/ \n")
                fn.write("cd rep"+str(i)+'.c'+c+" \n\n\n")

                fn.write('/usr/bin/time -v -o timeNmem.metilene_unsup python '+wd+'/metilene3/metilene3.py \
                        -i ../input/'+ab+'a10.metilene3.tsv \
                        -w 10 \
                        -wsup False \
                        -o ./metilene_unsup \
                        -t '+c+' -plot False \n')
                fn.write("cd ..;cp -rf rep"+str(i)+'.c'+c+" "+outdir+"/ \n")
                fn.write("cd rep"+str(i)+'.c'+c+" \n\n\n")

                
                fn.write("cd ..;cp -rf rep"+str(i)+'.c'+c+" "+outdir+"/\n\n\n")
    
    os.system("chmod +x "+mxqsh)
    os.system("mxqsub --processors=10 --runtime=7d --memory=500G --tmpdir=1000G \
                    --stdout="+mxqsh+".log \
                    --stderr="+mxqsh+".err \
                    "+mxqsh)

# %%
ab = '40_3_22_22'
run_all(ab)

# %%
ab = '15_5_10_10'
run_all(ab)


