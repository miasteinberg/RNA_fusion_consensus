# vim: ft=python
import os
import glob
import itertools
import pandas
from collections import defaultdict

configfile: 'fusion_config.yaml'
workdir: os.environ['PWD']
shell.executable('bash')

## This workflow takes a consensus approach with four fusion detection programs.
# Once candidate fusions have been identified, users can look closer at the output statistics to better identify real and false positives 
# (e.g. no spanning frags + no large anchor support = suspicious)


## ---- The parser may have to be customized for each run ---- ##
def parse_sampleID(filename):
    return filename.split('/')[-1].split('_')[0]

fastqs = glob.glob(config['fastq_dir'] + '/SC08063*fastq.gz')

d = defaultdict(list) # pair the R1 and R2 files by sample ID
for key, value in itertools.groupby(fastqs, parse_sampleID):
    d[key] += list(value)

sampleIDs = d.keys()

def input_files(wildcards):
    return sorted(d[wildcards.sampleID])

def sorted_fusion(gene1, gene2):
    return '|'.join(sorted([gene1, gene2]))


include: 'ericscript_Snakefile'
include: 'mapsplice_Snakefile'
include: 'starfusion_Snakefile'
include: 'chimerascan_Snakefile'


# These rules run on the host node and are not submitted to the cluster.
localrules: all


rule all:
    input:  
        'fusions_all_scores_merged.tsv'
 

# specifically coded for these four tools - will need modification to add more
# score the fusions for each tool and sample as the file is read and appended
# create one table for each tool with all samples in each table
rule build_tables:
    input: 
        eric = expand(config['es_calls'] + '/{sampleID}' + config['es_suff'], sampleID=sampleIDs),
        map = expand(config['ms_calls'] + '/{sampleID}' + config['ms_suff'], sampleID=sampleIDs),
        star = expand(config['sf_calls'] + '/{sampleID}' + config['sf_suff'], sampleID=sampleIDs),
        chim = expand(config['cs_calls'] + '/{sampleID}' + config['cs_suff'], sampleID=sampleIDs)
    output: 
        eric = 'tables/all_eric.txt',
        map = 'tables/all_map.txt',
        star = 'tables/all_star.txt',
        chim = 'tables/all_chim.txt'
    params: 
        span = 3, # this is how many supporting reads needed to be included
        mscols = config['ms_cols'] # list of column headers
    run:
        # create ES table
        es_dfs = []
        for fname in input.eric:
            df = pandas.read_table(fname)
            df['sampleID'] = fname.split('/')[-1].split('_')[0]
            df['ericscript'] = 1
            df['fusion'] = df.apply(lambda x: sorted_fusion(x.GeneName1, x.GeneName2), axis=1)
            df['supporting_reads'] = df.spanningreads
            df['total_support'] = df['supporting_reads'].groupby(df['fusion']).transform('sum') # combine fusions that are A|B and B|A
            df.drop_duplicates('fusion', inplace=True)  # only keep first row with fusion now that support reads are summed
            df = df[df.total_support >= params.span]    # remove fusions with too few support reads
            scores = list(range(1, len(df) + 1))
            scores.reverse()    # you want the fusions w/most reads getting highest score
            df.sort_values(by=['total_support'], ascending=[False], inplace=True)
            df['es_rank'] = scores
            df['es_score'] = df.es_rank.apply(lambda x: float(x)/len(df))   # percent scores for each fusion (1.0 is top fusion)
            es_dfs.append(df)
        es_df = pandas.concat(es_dfs)
        es_df.to_csv(output.eric, sep='\t')

        # create MS table
        ms_dfs = []
        for fname in input.map:
            df = pandas.read_table(fname, names=params.mscols, index_col=False)
            df['sampleID'] = fname.split('/')[-1].split('_')[0]
            df['mapsplice'] = 1
            df['gene1'] = df.annotated_gene_donor.apply(lambda x: x.split(',')[0])
            df['gene2'] = df.annotated_gene_acceptor.apply(lambda x: x.split(',')[0])
            df['fusion'] = df.apply(lambda x: sorted_fusion(x.gene1, x.gene2), axis=1)
            df['supporting_reads'] = df['coverage']
            df['total_support'] = df['supporting_reads'].groupby(df['fusion']).transform('sum') # combine fusions that are A|B and B|A
            df.drop_duplicates('fusion', inplace=True)  # only keep first row with fusion now that support reads are summed
            df = df[df.total_support >= params.span]    # remove fusions with too few support reads
            scores = list(range(1, len(df) + 1))
            scores.reverse()    # you want the fusions w/most reads getting highest score
            df.sort_values(by=['total_support'], ascending=[False], inplace=True)
            df['ms_rank'] = scores
            df['ms_score'] = df.ms_rank.apply(lambda x: float(x)/len(df))   # percent scores for each fusion (1.0 is top fusion)
            ms_dfs.append(df)
        ms_df = pandas.concat(ms_dfs)
        ms_df.to_csv(output.map, sep='\t')

        # create SF table
        sf_dfs = []
        for fname in input.star:
            df = pandas.read_table(fname)
            df['sampleID'] = fname.split('/')[-1].split('_')[0]
            df['starfusion'] = 1
            df['gene1'] = df['LeftGene'].apply(lambda x: x.split('^')[0])
            df['gene2'] = df['RightGene'].apply(lambda x: x.split('^')[0])
            df['fusion'] = df.apply(lambda x: sorted_fusion(x.gene1, x.gene2), axis=1)
            df['supporting_reads'] = df['SpanningFrags'] + df['JunctionReads']
            df['total_support'] = df['supporting_reads'].groupby(df['fusion']).transform('sum') # combine fusions that are A|B and B|A
            df.drop_duplicates('fusion', inplace=True)  # only keep first row with fusion now that support reads are summed
            df = df[df.total_support >= params.span]    # remove fusions with too few support reads
            scores = list(range(1, len(df) + 1))
            scores.reverse()    # you want the fusions w/most reads getting highest score
            df.sort_values(by=['total_support'], ascending=[False], inplace=True)
            df['sf_rank'] = scores
            df['sf_score'] = df.sf_rank.apply(lambda x: float(x)/len(df))   # percent scores for each fusion (1.0 is top fusion)
            sf_dfs.append(df)
        sf_df = pandas.concat(sf_dfs)
        sf_df.to_csv(output.star, sep='\t')

        # create CS table
        cs_dfs = []
        for fname in input.chim:
            df = pandas.read_table(fname)
            df['sampleID'] = fname.split('/')[-1].split('_')[0]
            df['chimerascan'] = 1
            df['fusion'] = df.apply(lambda x: sorted_fusion(x.genes5p, x.genes3p), axis=1)
            df['supporting_reads'] = df.spanning_frags
            df['total_support'] = df['supporting_reads'].groupby(df['fusion']).transform('sum') # combine fusions that are A|B and B|A
            df.drop_duplicates('fusion', inplace=True)  # only keep first row with fusion now that support reads are summed
            df = df[df.total_support >= params.span]    # remove fusions with too few support reads
            scores = list(range(1, len(df) + 1))
            scores.reverse()    # you want the fusions w/most reads getting highest score
            df.sort_values(by=['total_support'], ascending=[False], inplace=True)
            df['cs_rank'] = scores
            df['cs_score'] = df.cs_rank.apply(lambda x: float(x)/len(df))   # percent scores for each fusion (1.0 is top fusion)
            cs_dfs.append(df)
        cs_df = pandas.concat(cs_dfs)
        cs_df.to_csv(output.chim, sep='\t')


rule combine_scores:
    input: 
        eric = 'tables/all_eric.txt',
        map = 'tables/all_map.txt',
        star = 'tables/all_star.txt',
        chim = 'tables/all_chim.txt'

    output: 'tables/all_scores.txt'
    params: meta = config['metadata']
    run:
        meta = pandas.read_excel(params.meta, skiprows=2)
        meta.dropna(subset=['CGR Sample ID for seq'], inplace=True)
        meta['sampleID'] = meta['CGR Sample ID for seq'].apply(lambda x: x.rsplit('-', 1)[0])
        meta2 = meta[['TSS Subject ID', 'sampleID', 'Specimen Type ID']].copy()

        df_es = pandas.read_table(input.eric, sep='\t')
        df_ms = pandas.read_table(input.map, sep='\t')
        df_sf = pandas.read_table(input.star, sep='\t')
        df_cs = pandas.read_table(input.chim, sep='\t')

        cols = ['sampleID', 'fusion']
        df = df_es[cols + ['ericscript', 'es_score']].merge(df_ms[cols + ['mapsplice', 'ms_score']], on=cols, how='outer')
        df = df.merge(df_sf[cols + ['starfusion', 'sf_score']], on=cols, how='outer')
        df = df.merge(df_cs[cols + ['chimerascan', 'cs_score']], on=cols, how='outer')
        df.fillna(0, inplace=True)

        df['sum'] = df.ericscript + df.mapsplice + df.starfusion + df.chimerascan
        df['total_score'] = df.es_score + df.ms_score + df.sf_score + df.cs_score

        before = len(df)
        df = df.merge(meta2, on='sampleID', how='left') # add patient ID and tissue type
        after = len(df)
        df.to_csv(output[0], sep='\t', index=False)
        
        if before != after:
            shell('echo "Some samples were lost during the metadata merge (%d !+ %d). Check the metadata table for completeness." > LOST_SAMPLES.OHNOEZ' %(before, after))

        
rule merge_patient:
    input: rules.combine_scores.output
    output: 'fusions_all_scores_merged.tsv'
    run:        
        df = pandas.read_table(input[0], sep='\t')
        dt = df[df['Specimen Type ID'] == 'Primary Tumor'].copy()
        dn = df[df['Specimen Type ID'] == 'Solid Tissue Normal'].copy()
        dd = dt.merge(dn, how='outer', on=['TSS Subject ID', 'fusion'], suffixes=('_t', '_n'))
        ddd = dd[['TSS Subject ID', 'fusion', 'sum_n', 'total_score_n', 'sum_t', 'total_score_t']].sort_values(['sum_t', 'total_score_t', 'sum_n', 'total_score_n'], ascending=False)
        ddd.to_csv(output[0], sep='\t', index=False)


