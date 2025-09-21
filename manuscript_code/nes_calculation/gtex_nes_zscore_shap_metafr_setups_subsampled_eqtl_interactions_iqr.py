import random
import numpy as np
import pandas as pd
import os
from pathlib import Path
#from pybedtools import BedTool
import pybedtools
import gseapy as gp
import sys
from collections import defaultdict
from itertools import chain
import gzip
sys.path.append("/home/dhecker/FuFis/src/")
import GTEx_eQTLReader
import GTF_Processing
from scipy.stats import median_abs_deviation
from sklearn.preprocessing import RobustScaler
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
import seaborn as sns

#### comparison of different metafr setups: zscore SHAP #####

bedtools_tempdir = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/pybedtools_tmp"
os.makedirs(bedtools_tempdir, exist_ok=True)
pybedtools.set_tempdir(bedtools_tempdir)

from pybedtools import BedTool

#random.seed(1234)

base_folder = "/projects/apog/work/IHEC/ValidateInteractions/eQTL_Validation/"
hg38_annotation = "/projects/abcpp/work/base_data/gencode.v38.annotation.gtf"
chrom_sizes = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/chrom_sizes_BSgenome_Hsapiens_UCSC_hg38.txt"
#TODO: use only genes for which all setups could be trained!!

#interactions_sc = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions"
#interactions_mc100 = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions"
#interactions_mc50 = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions"
#interactions_mc150 = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions"
#interactions_scarlink = "/projects/single_cell_stitchit/work/scarlink/pbmc_output_100pb_bins_train/gene_linked_tiles_celltypes.csv"

#nes_out_path_sc = '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood.tsv'
#nes_out_path_mc100 = '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood.tsv'
#nes_out_path_mc50 = '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc50_nes_whole_blood.tsv'
#nes_out_path_mc150 = '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc150_nes_whole_blood_extended_bins_200bp_top100000.tsv'

#comparison metafr setups: zscore 
""" config = {
    
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_abs_zscore_distance.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc50_nes_whole_blood_ext200bp_top100000_abs_zscore_distance.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood_ext200bp_top100000_abs_zscore_distance.tsv'
    }, 
    'mc150': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc150_nes_whole_blood_ext200bp_top100000_abs_zscore_distance.tsv'
    }
}

#comparison metafr setups: shap
config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_top100000_abs_shap.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc50_nes_whole_blood_top100000_abs_shap.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood_top100000_abs_shap.tsv'
    }, 
    'mc150': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc150_nes_whole_blood_top100000_abs_shap.tsv'
    }
} """

#500bp

#shap
"""
config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/sc_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc50_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc100_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    }, 
    'mc150': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc150_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    }
}"""

config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_variable_1000_genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_variable_1000_genes/eqtl_formatted_interactions/shap/sc_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_1000_variable_genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_1000_variable_genes/eqtl_formatted_interactions/shap/mc50_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_1000_variable_genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_1000_variable_genes/eqtl_formatted_interactions/shap/mc100_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    }, 
    'mc150': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_variable_1000_genes/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_variable_1000_genes/eqtl_formatted_interactions/shap/mc150_nes_whole_blood_top100000_abs_shap_distance_iqr.tsv'
    }
}

""" #zscore
config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/zscore",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/sc_nes_whole_blood_top100000_abs_shap_zscore.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/zscore",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc50_nes_whole_blood_top100000_abs_shap_zscore.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/zscore",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc100_nes_whole_blood_top100000_abs_shap_zscore.tsv'
    }, 
    'mc150': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/zscore",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc150_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions/shap/mc150_nes_whole_blood_top100000_abs_shap_zscore.tsv'
    }
}


#comparison zscore vs shap for single-cell metafr setup
config = {
    'sc_zscore': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_abs_zscore_shap.tsv'
    },
    'sc_shap': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_abs_shap.tsv'
    }
}
 """

'''
config = {
    'mc100_100bp': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood_100bp_top100000.tsv'
    },
    'mc100_250bp': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB250bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB250bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood_250bp_top100000.tsv'
    },
    'mc100_500bp': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB500bp_binned_lsi_activity_q0_no_threshold_global/mc100_nes_whole_blood_500bp_top100000.tsv'
    }
}
'''

'''
config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes/sc_nes_whole_blood_ext200bp_top100000_top2000genes.tsv'
    },
    'mc50': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc50_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/mc50_nes_whole_blood_ext200bp_top100000_top2000genes.tsv'
    },
    'mc100': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/mc100_1MB100bp_binned_lsi_activity_q0_no_threshold_global_top2000genes/mc100_nes_whole_blood_ext200bp_top100000_top2000genes.tsv'
    }
}
'''
'''
config = {
    'sc': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/eqtl_formatted_interactions",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_no_abs_ranking.tsv'
    }
}

config = {
    'sc_zscore': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_ct_unspecific/eqtl_formatted_interactions/zscore",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_zscore_ct_unspecific.tsv'
    },
    'sc_shap': {
        'interactions_path': "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_ct_unspecific/eqtl_formatted_interactions/shap",
        'output_path': '/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global/sc_nes_whole_blood_ext200bp_top100000_shap_ct_unspecific.tsv'
    }
}
'''

n_top = 100000

#Interaction                            Score        
#ENSG00000179055\t32745308\t32755308    0.740747
bin_size=100
ext_factor = 2
extension = ext_factor*bin_size

#TODO: extension as flag
extension_flag = False #TODO: 500 bp bins -> extension?

#TODO: only use common genes between setups
def prepare_interactions(interactions_file):
    print("load interactions...")
    celltype_dfs = {}
    for file in Path(interactions_file).glob("*_interactions_all_genes.tsv"):
        cell_type = file.stem.split('_interactions')[0].replace('zscore_shap_', '').replace('shap_', '')
        cell_type = cell_type.replace(' ', '_').strip('_')
        df = pd.read_csv(file, sep='\t')
        celltype_dfs[cell_type] = df 
    return(celltype_dfs)

def subsample_hits(eqtls, target_size, setup_name, method_name):
    print("subsample interactions...")
    seed = hash(f"{setup_name}_{method_name}") % (2**32) 
    random.seed(seed)
    return BedTool(random.sample(list(eqtls), target_size))

def calculate_top_interactions(df):
    print("calculate top interactions...")
    top_df = (
            df.copy()
            .sort_values(by=['score'], ascending=False, ignore_index=False, key=abs) #TODO: try sorting positive and negative interactions
            .head(n_top)
        )
    return(top_df)

#TODO: include distance metric
def calculate_top_interactions_distance_weighted(df, n_top, dist_df):
    print("calculate top interactions...")
    weights = np.exp(-dist_df['distance_to_tss'] / 200000) #200kb fixed from GEEES and scReg
    weighted_scores = np.abs(df['score']) * weights
    top_df = (
        df.copy()
        .assign(weighted_score=weighted_scores)
        .sort_values('weighted_score', ascending=False)
        .head(n_top)
    )
    return(top_df)

#TODO: use robust z-score and include distance metric
def calculate_top_interactions_distance_robust_zscore_weighted(df, n_top, dist_df, epsilon=1e-6, min_mad=1e-6):
    print("calculate top interactions...")
    def compute_mad_zscore(shap_values):
        shap_regularized = shap_values + epsilon  # Avoid MAD=0
        med = np.median(shap_regularized)
        mad = median_abs_deviation(shap_regularized, scale='normal')  # 1.4826 scaling
        mad = max(mad, min_mad)
        return (shap_values - med) / mad  
    if 'gene_id' not in df.columns:
        merged_df = pd.merge(
            df, 
            dist_df[['interaction', 'gene_id', 'distance_to_tss']], 
            on='interaction',
            how='left' #same interactions in both dfs
        )
        result_df = merged_df 
    else: result_df = df.copy() #first entry in df dictionary       
    result_df['robust_z'] = result_df.groupby('gene_id', group_keys=False)['score'].apply(compute_mad_zscore)
    weights = np.exp(-result_df['distance_to_tss'] / 200000) #200kb fixed from GEEES and scReg
    result_df['weighted_score'] = np.abs(result_df['robust_z']) * weights
    top_df = result_df.sort_values('weighted_score', ascending=False).head(n_top).copy()
    return(top_df)

def calculate_top_interactions_distance_robust_iqr_scaling_weighted(df, n_top, dist_df):
    print("calculate top interactions...")
    scaler = RobustScaler(quantile_range=(25, 75))
    if 'gene_id' not in df.columns:
        merged_df = pd.merge(
            df, 
            dist_df[['interaction', 'gene_id', 'distance_to_tss']], 
            on='interaction',
            how='left' #same interactions in both dfs
        )
        result_df = merged_df 
    else: result_df = df.copy() #first entry in df dictionary       
    def scale_iqr(gene):
        gene['robust_z'] = scaler.fit_transform(gene[['score']]).flatten()
        return gene
    result_df = result_df.groupby('gene_id', group_keys=False).apply(scale_iqr)
    weights = np.exp(-result_df['distance_to_tss'] / 200000) #200kb fixed from GEEES and scReg
    result_df['weighted_score'] = np.abs(result_df['robust_z']) * weights
    top_df = result_df.sort_values('weighted_score', ascending=False).head(n_top).copy()
    return(top_df)

def formatting_extended_interactions(df, extended_regions_dict, full_df):
    print("calculate mean scores across overlapping interactions...")
    results = []
    full_score_series = full_df['score']
    for interaction in df['interaction']:
        gene, start, end = interaction.split('\t')
        ext_data = extended_regions_dict.get((gene,start,end))
        extended_region = ext_data['extended_region']
        #print(extended_region)
        overlapping_indices = ext_data['overlapping_bins']
        #print(overlapping_indices)
        #overlapping_scores = score_series.iloc[overlapping_indices]
        overlapping_scores = full_score_series.loc[overlapping_indices]
        #print(overlapping_scores)
        if len(overlapping_scores) > 0:
            results.append({
                'interaction': f"{gene}\t{extended_region[1]}\t{extended_region[2]}",
                'score': overlapping_scores.abs().mean() #TODO: keep sign? different signs overlapping regions? check
            })
    mean_score_df = pd.DataFrame(results)
    ext_interactions_bed = BedTool.from_dataframe(
        mean_score_df['interaction'].str.split('\t', expand=True)
        #.rename(columns={0: 'gene', 1: 'start', 2: 'end'})
    )
    return mean_score_df, ext_interactions_bed

'''
def formatting_extended_interactions(df, extension_cache):
    extended_interactions = []
    for interaction in df['interaction']:
            gene = interaction.split('\t')[0]
            #print(gene)
            start = interaction.split('\t')[1]
            #print(start)
            end = interaction.split('\t')[2]
            #print(end)
            ext_interaction = extension_cache.get((gene,start,end))
            #print(ext_interaction)
            #extended_interactions.append(ext_interaction)
            ext_interaction_formatted = (
                    gene + '\t' + 
                    str(ext_interaction[1]) + '\t' + 
                    str(ext_interaction[2])
            )
            extended_interactions.append(ext_interaction_formatted)       
    ext_interactions_bed = BedTool('\n'.join(pd.Series(extended_interactions)), from_string=True) 
    return([extended_interactions,ext_interactions_bed])           
'''

def intersect_eqtls(interactions_bed, method):
    print("intersect eqtls...")
    eqtl_hits_ct = interactions_bed.intersect(gtex_eqtl_beds['Whole_Blood'][method], u=True) #positive interactions
    eqtls_ct = gtex_eqtl_beds['Whole_Blood'][method].intersect(interactions_bed, u=True) #positive eqtls
    return(eqtl_hits_ct, eqtls_ct)

#gene_body_bed(gtf_file, gene_set=set(), dict_only=True)
def get_genes(df):
    genes = df.iloc[:, 0].str.split('\t').str[0]
    unique_genes = set(genes)
    return(unique_genes)

'''
def extend_interactions(df, chr_sizes, g_coords):
    extended_regions_dict = {}
    for interaction in df['interaction']:
        gene = interaction.split('\t')[0]
        start = interaction.split('\t')[1]
        end = interaction.split('\t')[2]
        chrom = g_coords[gene][0]
        chrom_size = chr_sizes[chrom]
        new_start = max(1, int(start) - extension)
        new_end = min(chrom_size, int(end) + extension)
        extended_regions_dict[(gene, str(start), str(end))] = ((gene, new_start, new_end))
    return(extended_regions_dict)
'''
def extend_interactions(df, chr_sizes, g_coords, extension, ext_factor):
    print("extending interactions...")
    extended_regions_dict = {}
    bin_indices = df.index.tolist() #TODO: max index?
    max_idx = len(bin_indices) - 1
    for idx, row in df.iterrows():
        interaction = row['interaction']
        gene = interaction.split('\t')[0]
        start = interaction.split('\t')[1]
        end = interaction.split('\t')[2]
        chrom = g_coords[gene][0]
        chrom_size = chr_sizes[chrom]
        new_start = max(1, int(start) - extension)
        #print(new_start)
        new_end = min(chrom_size, int(end) + extension)
        #print(new_end)
        start_idx = max(0, int(idx) - ext_factor) #row index
        end_idx = min(max_idx, int(idx) + ext_factor)
        overlapping_idx = list(range(int(start_idx), int(end_idx) + 1)) 
        #print(overlapping_idx)       
        extended_regions_dict[(gene, start, end)] = {
            'extended_region': (gene, new_start, new_end),
            'overlapping_bins': overlapping_idx
        }
    return(extended_regions_dict)
    
#all shap dfs contain same interactions -> use random df for storing the extended region coordinates
def extension_caching(df, chrom_sizes, extension, ext_factor):
    print("create extension cache for all interactions...")
    #ex_df = celltype_dfs['CD14_Mono']
    uniq_genes = get_genes(df)
    gene_coords = GTF_Processing.gene_body_bed(hg38_annotation, gene_set=uniq_genes, dict_only=True)
    #chrom_sizes_df = pd.read_csv(chrom_sizes, sep='\t', header=0, names=['chrom', 'size'])
    #chrom_sizes = dict(zip(chrom_sizes_df['chrom'], chrom_sizes_df['size']))
    extension_cache = extend_interactions(df, chrom_sizes, gene_coords, extension, ext_factor)
    return(extension_cache)

 
def calculate_nes(ext_interactions_df, eqtl_hits_ct): #celltype_df, #overlapping interactions
    print("calculate NES...")
    if len(eqtl_hits_ct) > 0:
                re_res = gp.prerank(rnk=ext_interactions_df, # CARE they sort themselves again by the first column.
                        gene_sets={'Interactions with overlap': set(['\t'.join(x.fields[:3]) for x in eqtl_hits_ct])},
                        min_size=1,
                        max_size=len(gtex_eqtl_beds['Whole_Blood']['DAP-G'])+1,  # Max hits. Doesn't affect the output though.
                        permutation_num=100, # Reduce number to speed up testing.
                        outdir=None,
                        seed=1234,
                        threads=1,
                        weight=0,
                        no_plot=True,
                        verbose=False,
                        ).res2d       
    return(re_res)

## Dennis function ##
def gene_window_bed(gtf_file, extend=200, gene_set=set(), tss_type='5', dict_only=False, merge=False,
                    open_regions=False):
    """
    Based on a gtf file fetches all or the most 5' TSS for all genes, and returns a BedTool object with windows
    around the TSS, expanding by 'extend' in each direction, resulting in a total window size of 2*'extend'+1.
    Alternatively gives a dictionary with the TSS, also containing the number of transcripts, gene name and strand.
    The BedTools intervals will be 0-based, the TSS in the dictionary still 1-based like in the gtf-file.
    Care: removes the .-suffixes from all gene IDs.

    Args:
        gtf_file: gtf-file in GENCODE's format, can be gzipped.
        extend: Number of base pairs to extend the TSS in each direction. 200 means a window of size 401.
        gene_set: Set of Ensembl IDs or gene names or mix of both to limit the output to. If empty, return for all
            genes in the annotation.
        tss_type: "5" to get only the 5' TSS or "all" to get all unique TSS of all transcripts in the gtf-file.
        dict_only: Returns a dictionary instead of a BedTool's object.
        merge: If True, merges all overlapping promoters of the same gene into one row in the BedTool's object.
        open_regions: Optional bed file or BedTools' object, only overlapping parts of promoters will be kept for the
            BedTool's object. Can therefore be used to find genes whose promoter overlap a set of peaks, for example to
            find genes that are accessible.
    """
    if tss_type == '5':
        identifier = 'gene'
    elif tss_type == 'all':
        identifier = 'transcript'
    if gtf_file.endswith('.gz'):
        file_opener = gzip.open(gtf_file, 'rt')
    else:
        file_opener = open(gtf_file)
    if gene_set is not None:
        gene_set = set([g.split('.')[0] for g in gene_set])
    tss_locs = {}
    with file_opener as gtf_in:
        for entry in gtf_in:
            if not entry.startswith('#') and entry.split('\t')[2] == identifier:
                line = entry.strip().split('\t')
                # Some gene IDs are non-unique if they have a _PAR_Y version.
                if not line[8].split('gene_id "')[-1].split('";')[0].endswith("_PAR_Y"):
                    this_gene = line[8].split('gene_id "')[-1].split('";')[0].split('.')[0]
                    gene_name = line[8].split('gene_name "')[-1].split('";')[0]
                    if not gene_set or this_gene in gene_set or gene_name in gene_set:
                        if this_gene not in tss_locs:
                            tss_locs[this_gene] = {'chr': None, 'tss': set(), '#transcripts': 0}
                        tss_locs[this_gene]['chr'] = line[0]
                        tss_locs[this_gene]['name'] = gene_name
                        if line[6] == '+':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] > int(line[3])):
                                tss_locs[this_gene]['tss'] = {int(line[3])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[3]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '+'
                        if line[6] == '-':
                            if identifier == 'gene' and (not tss_locs[this_gene]['tss'] or list(tss_locs[this_gene]['tss'])[0] < int(line[4])):
                                tss_locs[this_gene]['tss'] = {int(line[4])}
                            elif identifier == 'transcript':
                                tss_locs[this_gene]['tss'].add(int(line[4]))
                                tss_locs[this_gene]['#transcripts'] += 1
                            tss_locs[this_gene]['strand'] = '-'
    if dict_only:
        return tss_locs
    promoter_bed = BedTool('\n'.join(chain(*[[vals['chr'] + '\t' + str(max([0, tss - int(extend) - 1])) + '\t' +
                                             str(tss + int(extend)) + '\t' + g + '\t.\t' + vals['strand'] for tss in vals['tss']]
                                             for g, vals in tss_locs.items()])), from_string=True)
    if open_regions and str(open_regions).lower() != "false":
        if type(open_regions) == str:
            open_regions = BedTool('\n'.join(['\t'.join(x.strip().split('\t')[:3]) for x
                                              in open(open_regions).readlines() if not x.startswith('#')]), from_string=True)
        promoter_bed = promoter_bed.intersect(open_regions)
    if merge:  # Flip the chr and Ensembl ID column to merge promoter of the same gene, and afterwards flip again.
        # Due to using the Ensembl ID as the first coordinate we don't need to consider the strand for merging (relying
        # on the annotations putting a gene only on one strand exclusively).
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True).sort().merge(c=[4, 5, 6], o='distinct')
        promoter_bed = BedTool('\n'.join(['\t'.join([x.fields[3], x.fields[1], x.fields[2], x.fields[0], x.fields[4], x.fields[5]]) for x in promoter_bed]), from_string=True)
    return promoter_bed

#results_df.to_csv(nes_out_path, sep='\t', index=False) 
#if __name__ == "__main__":
gtex_eqtl_beds, gtex_eqtl_genes = GTEx_eQTLReader.get_eqtls(gtex_folder=base_folder, gtex_tissues=False, 
                                                        hg38_annotation=hg38_annotation, max_distance=1000000)  

chrom_sizes_df = pd.read_csv(chrom_sizes, sep='\t', header=0, names=['chrom', 'size'])
chrom_sizes = dict(zip(chrom_sizes_df['chrom'], chrom_sizes_df['size']))

eqtl_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
min_counts = defaultdict(lambda: defaultdict(int)) 

#find minimum number of overlapping eqtls across methods
#min_eqtl_hits = {}
common_celltypes = None
gtex_methods = ['DAP-G']#, 'CaVEMaN', 'CAVIAR']

ex_df = {}
ct_df_dict = {}

if extension_flag:
    ext_interactions_cache = {}

for setup, files in config.items():
    ct_ext_cache = {}
    interactions = files.get('interactions_path')
    print(interactions)
    ct_df_dict[setup] = prepare_interactions(interactions) 
    setup_celltypes = set(ct_df_dict[setup].keys())
    if common_celltypes is None:
        common_celltypes = setup_celltypes
    else:
        common_celltypes.intersection_update(setup_celltypes)
    ct = next(iter(ct_df_dict[setup]))
    ct_df = ct_df_dict[setup].get(ct)     
    if extension_flag:
        #ct = next(iter(ct_df_dict))
        #ct_df = ct_df_dict.get(ct) #chose random celltype to calculate bin interactions, TODO: optimize -> equal for all mc setups
        ct_ext_cache = extension_caching(ct_df, chrom_sizes, extension, ext_factor) # extend all interactions
    #TODO: genes, regions equal for all cell types: calculate distances between 5' tss and genomic regions once
    coord_arr = ct_df.iloc[:, 0].to_numpy() #actually the same for all setups -> #TODO: just calculate once
    split_coords = np.array([x.split('\t') for x in coord_arr])
    gene_ids = split_coords[:, 0] 
    starts = split_coords[:, 1].astype(int)  
    ends = split_coords[:, 2].astype(int)  
    target_genes = np.unique(gene_ids)
    #print(target_genes)
    tss_coords = gene_window_bed(gtf_file=hg38_annotation, extend=0, gene_set=target_genes, tss_type='5', dict_only=True)
    print(tss_coords) 
    tss_positions = np.array([next(iter(tss_coords[gene]['tss'])) for gene in gene_ids], dtype=np.int32)
    dist_start = np.abs(starts - tss_positions)
    dist_end = np.abs(ends - tss_positions)
    ct_df['distance_to_tss'] = np.minimum(dist_start, dist_end)
    ct_df['gene_id'] = gene_ids
    ex_df[setup] = ct_df
    for cell_type, df in ct_df_dict[setup].items():
        print(cell_type)
        #print(df)
        if cell_type not in common_celltypes: #only use common cell types
            print("celltype not in common celltypes - skip")
            continue
        #top_df = calculate_top_interactions(df) #top interactions = absolute shap value, cut-off on #top interactions
        #top_df =  calculate_top_interactions_distance_weighted(df, n_top, ex_df[setup])
        #top_df = calculate_top_interactions_distance_robust_zscore_weighted(df, n_top, ex_df[setup], epsilon=1e-6, min_mad=1e-6) #use robust z-score
        top_df = calculate_top_interactions_distance_robust_iqr_scaling_weighted(df, n_top, ex_df[setup]) #use robust z-score
        if extension_flag:
            mean_score_df, top_interactions_bed = formatting_extended_interactions(top_df, ct_ext_cache, df) # return([extended_interactions,ext_interactions_bed])
            ext_interactions_cache[(setup, cell_type)] = (mean_score_df, top_interactions_bed)     
        else:
            top_interactions_bed = BedTool('\n'.join(top_df['interaction'].str.split('\t').apply(lambda x: '\t'.join(x[:3]))), from_string=True) #use top interactions without extension
        method_eqtl_hits = {}
        for method in gtex_methods:
            eqtl_hits, eqtls = intersect_eqtls(top_interactions_bed, method) #interactions, eqtls (positive)
            eqtl_counts[setup][cell_type][method] = len(eqtls)

all_cell_types = set() 
for setup in eqtl_counts:
    all_cell_types.update(eqtl_counts[setup].keys())
for cell_type in all_cell_types:
    for method in gtex_methods:
        counts = []
        for setup in eqtl_counts:
            if cell_type in eqtl_counts[setup]:
                counts.append(eqtl_counts[setup][cell_type][method])
        if counts:
            min_counts[cell_type][method] = min(counts)    

for setup, files in config.items():
    interactions = files.get('interactions_path')
    print(interactions)
    out_path = files.get('output_path')
    print(out_path)
    results = []
    for cell_type, df in ct_df_dict[setup].items():
        print(cell_type)
        print(df)
        if cell_type not in common_celltypes: #only use common cell types
            print("celltype not in common celltypes - skip")
            continue
        #top_df = calculate_top_interactions(df)
        #top_df =  calculate_top_interactions_distance_robust_zscore_weighted(df, n_top, ex_df[setup], epsilon=1e-6, min_mad=1e-6) #use robust z-score
        top_df = calculate_top_interactions_distance_robust_iqr_scaling_weighted(df, n_top, ex_df[setup]) #use robust z-score
        try:
            fig, ax = plt.subplots(figsize=(8,5))
            sns.kdeplot(top_df["robust_z"], bw_adjust=0.5, fill=True, color='royalblue', alpha=0.6, linewidth=2, ax=ax)
            ax.set_title("Distribution of Robust-scaled SHAP Values", fontsize=16, weight='bold')
            ax.set_xlabel("Robust-scaled SHAP", fontsize=14)
            ax.set_ylabel("Density", fontsize=14)
            ax.grid(True, linestyle='--', alpha=0.4)
            sns.despine(ax=ax)
            plt.tight_layout()
            save_path = os.path.join(interactions, f"{cell_type}_robust_shap_density_plot.pdf")
            fig.savefig(save_path)
            plt.close(fig)
            print(f" - Saved plot to {save_path}")
        except Exception as e:
            print(f" - Failed to plot for {cell_type}: {e}")
        """"
        os.makedirs(interactions, exist_ok=True) #obsolete
        ### plot distribution of robust-scaled shap values ###
        plt.figure(figsize=(8,5))
        sns.kdeplot(top_df["robust_z"], bw_adjust=0.5, fill=True, color='royalblue', alpha=0.6, linewidth=2)

        plt.title("Distribution of Robust-scaled SHAP Values", fontsize=16, weight='bold')
        plt.xlabel("Robust-scaled SHAP", fontsize=14)
        plt.ylabel("Density", fontsize=14)
        plt.grid(True, linestyle='--', alpha=0.4)
        sns.despine()  # clean top/right borders

        plt.tight_layout()
        plt.savefig(os.path.join(interactions, f"{cell_type}_robust_shap_density_plot.pdf"))
        print(plt.savefig(os.path.join(interactions, f"{cell_type}_robust_shap_density_plot.pdf")))
        plt.close()
        """
        ##########
        #top_df = top_df.drop('score', axis=1)
        print(top_df)
        if extension_flag:
            mean_score_df, top_interactions_bed =  ext_interactions_cache[(setup, cell_type)]# return([extended_interactions,ext_interactions_bed])   
            #ext_dict = {'interaction' : top_interactions}
            #ext_top_df = pd.DataFrame(ext_dict)
            #ext_top_df['absolute_score'] = top_df['absolute_score']
        else:
            top_interactions_bed = BedTool('\n'.join(top_df['interaction'].str.split('\t').apply(lambda x: '\t'.join(x[:3]))), from_string=True)         
        for method in gtex_methods:
            eqtl_hits, eqtls = intersect_eqtls(top_interactions_bed, method)
            eqtls_subsampled = subsample_hits(eqtls, target_size=min_counts[cell_type][method], setup_name=setup, method_name=method)    
            eqtls_subsampled_bed = BedTool.from_dataframe(eqtls_subsampled.to_dataframe())
            interactions_with_subsampled_eqtls = top_interactions_bed.intersect(eqtls_subsampled_bed, u=True)
            #eqtls_subsampled = subsample_hits(eqtls, min_counts[cell_type][method])
            if extension_flag:
                ges_result = calculate_nes(mean_score_df, interactions_with_subsampled_eqtls)
            else:
                ges_result = calculate_nes(top_df, interactions_with_subsampled_eqtls)
            ges_result['CellType'] = cell_type
            ges_result['Method'] = method
            ges_result['Num_eQTL_hits'] = len(interactions_with_subsampled_eqtls)
            print(ges_result)
            results.append(ges_result)
    results_df = pd.concat(results)
    cols = ['CellType', 'Method', 'Num_eQTL_hits', 'ES', 'NES', 
    'NOM p-val', 'FDR q-val', 'FWER p-val', 'Tag %', 
    'Gene %', 'Lead_genes', 'Name']
    results_df = results_df[cols]
    print(results_df)
    results_df.to_csv(out_path, sep='\t', index=False)

        



"""
if __name__ == "__main__":
    gtex_eqtl_beds, gtex_eqtl_genes = GTEx_eQTLReader.get_eqtls(gtex_folder=base_folder, gtex_tissues=False, 
                                                            hg38_annotation=hg38_annotation, max_distance=1000000)  
    
    chrom_sizes_df = pd.read_csv(chrom_sizes, sep='\t', header=0, names=['chrom', 'size'])
    chrom_sizes = dict(zip(chrom_sizes_df['chrom'], chrom_sizes_df['size']))

    eqtl_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    min_counts = defaultdict(lambda: defaultdict(int)) 

    #find minimum number of overlapping eqtls across methods
    #min_eqtl_hits = {}
    common_celltypes = None
    gtex_methods = ['DAP-G', 'CaVEMaN', 'CAVIAR']

    ex_df = []

    if extension_flag:
        ext_interactions_cache = {}

    for setup, files in config.items():
        ct_ext_cache = {}
        interactions = files.get('interactions_path')
        print(interactions)
        ct_df_dict = prepare_interactions(interactions) 
        setup_celltypes = set(ct_df_dict.keys())
        if common_celltypes is None:
            common_celltypes = setup_celltypes
        else:
            common_celltypes.intersection_update(setup_celltypes)
        ct = next(iter(ct_df_dict))
        ct_df = ct_df_dict.get(ct)     
        if extension_flag:
            #ct = next(iter(ct_df_dict))
            #ct_df = ct_df_dict.get(ct) #chose random celltype to calculate bin interactions, TODO: optimize -> equal for all mc setups
            ct_ext_cache = extension_caching(ct_df, chrom_sizes, extension, ext_factor) # extend all interactions
        #TODO: genes, regions equal for all cell types: calculate distances between 5' tss and genomic regions once
        coord_arr = ct_df.iloc[:, 0].to_numpy() #actually the same for all setups -> #TODO: just calculate once
        split_coords = np.array([x.split('\t') for x in coord_arr])
        gene_ids = split_coords[:, 0] 
        starts = split_coords[:, 1].astype(int)  
        ends = split_coords[:, 2].astype(int)  
        target_genes = np.unique(gene_ids)
        #print(target_genes)
        tss_coords = gene_window_bed(gtf_file=hg38_annotation, extend=0, gene_set=target_genes, tss_type='5', dict_only=True)
        print(tss_coords) 
        tss_positions = np.array([next(iter(tss_coords[gene]['tss'])) for gene in gene_ids], dtype=np.int32)
        dist_start = np.abs(starts - tss_positions)
        dist_end = np.abs(ends - tss_positions)
        ct_df['distance_to_tss'] = np.minimum(dist_start, dist_end)
        ex_df = ct_df
        for cell_type, df in ct_df_dict.items():
            print(cell_type)
            #print(df)
            if cell_type not in common_celltypes: #only use common cell types
                print("celltype not in common celltypes - skip")
                continue
            #top_df = calculate_top_interactions(df) #top interactions = absolute shap value, cut-off on #top interactions
            top_df =  calculate_top_interactions_distance_weighted(df, n_top, ct_df)
            if extension_flag:
                mean_score_df, top_interactions_bed = formatting_extended_interactions(top_df, ct_ext_cache, df) # return([extended_interactions,ext_interactions_bed])
                ext_interactions_cache[(setup, cell_type)] = (mean_score_df, top_interactions_bed)     
            else:
                top_interactions_bed = BedTool('\n'.join(top_df['interaction'].str.split('\t').apply(lambda x: '\t'.join(x[:3]))), from_string=True) #use top interactions without extension
            method_eqtl_hits = {}
            for method in gtex_methods:
                eqtl_hits, eqtls = intersect_eqtls(top_interactions_bed, method) #interactions, eqtls (positive)
                eqtl_counts[setup][cell_type][method] = len(eqtls)

    all_cell_types = set() 
    for setup in eqtl_counts:
        all_cell_types.update(eqtl_counts[setup].keys())
    for cell_type in all_cell_types:
        for method in gtex_methods:
            counts = []
            for setup in eqtl_counts:
                if cell_type in eqtl_counts[setup]:
                    counts.append(eqtl_counts[setup][cell_type][method])
            if counts:
                min_counts[cell_type][method] = min(counts)    

    for setup, files in config.items():
        interactions = files.get('interactions_path')
        print(interactions)
        out_path = files.get('output_path')
        print(out_path)
        results = []
        for cell_type, df in ct_df_dict.items():
            print(cell_type)
            print(df)
            if cell_type not in common_celltypes: #only use common cell types
                print("celltype not in common celltypes - skip")
                continue
            #top_df = calculate_top_interactions(df)
            top_df =  calculate_top_interactions_distance_weighted(df, n_top, ex_df)
            #top_df = top_df.drop('score', axis=1)
            print(top_df)
            if extension_flag:
                mean_score_df, top_interactions_bed =  ext_interactions_cache[(setup, cell_type)]# return([extended_interactions,ext_interactions_bed])   
                #ext_dict = {'interaction' : top_interactions}
                #ext_top_df = pd.DataFrame(ext_dict)
                #ext_top_df['absolute_score'] = top_df['absolute_score']
            else:
                top_interactions_bed = BedTool('\n'.join(top_df['interaction'].str.split('\t').apply(lambda x: '\t'.join(x[:3]))), from_string=True)         
            for method in gtex_methods:
                eqtl_hits, eqtls = intersect_eqtls(top_interactions_bed, method)
                eqtls_subsampled = subsample_hits(eqtls, target_size=min_counts[cell_type][method], setup_name=setup, method_name=method)    
                eqtls_subsampled_bed = BedTool.from_dataframe(eqtls_subsampled.to_dataframe())
                interactions_with_subsampled_eqtls = top_interactions_bed.intersect(eqtls_subsampled_bed, u=True)
                #eqtls_subsampled = subsample_hits(eqtls, min_counts[cell_type][method])
                if extension_flag:
                    ges_result = calculate_nes(mean_score_df, interactions_with_subsampled_eqtls)
                else:
                    ges_result = calculate_nes(top_df, interactions_with_subsampled_eqtls)
                ges_result['CellType'] = cell_type
                ges_result['Method'] = method
                ges_result['Num_eQTL_hits'] = len(interactions_with_subsampled_eqtls)
                print(ges_result)
                results.append(ges_result)
        results_df = pd.concat(results)
        cols = ['CellType', 'Method', 'Num_eQTL_hits', 'ES', 'NES', 
        'NOM p-val', 'FDR q-val', 'FWER p-val', 'Tag %', 
        'Gene %', 'Lead_genes']
        results_df = results_df[cols]
        print(results_df)
        results_df.to_csv(out_path, sep='\t', index=False)
        

"""