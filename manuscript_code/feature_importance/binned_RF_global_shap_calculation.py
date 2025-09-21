import os
import pandas as pd
import joblib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import zscore
import shap
import os
from joblib import Parallel, delayed

#marginal global SHAP values

os.environ["OMP_NUM_THREADS"] = "10"
os.environ["MKL_NUM_THREADS"] = "10"

n_jobs=20

#output_path = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top2000genes"
output_path = "/projects/single_cell_stitchit/work/binned_RF/shap_analysis/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_global_top_variable_1000_genes"
os.makedirs(output_path, exist_ok=True)

log_file = f"{output_path}/shap_error_log_global.txt"

# no feature selection: all possible regions for comparing to scarlink shap
#model_path = "/projects/single_cell_stitchit/work/binned_RF/output/10x_pbmc_10k/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling/models"
model_path = "/projects/single_cell_stitchit/work/binned_RF/output/10x_pbmc_10k/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling/models"
#fmat_path = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/single_cell_1MB100bp_binned_lsi_activity_q0"
fmat_path = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/single_cell_1MB500bp_binned_lsi_activity_q0_test"
#performance_file = "/projects/single_cell_stitchit/work/binned_RF/output/10x_pbmc_10k/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling/results/sc_1MB100bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling_performance_correlation.txt"
#performance_file = "/projects/single_cell_stitchit/work/binned_RF/output/10x_pbmc_10k/sc_1MB500bp_binned_lsi_activity_q100_no_threshold_test_partition_no_scaling/sc_1MB500bp_binned_lsi_activity_q0_no_threshold_performance_correlation.txt"
gene_file = "/projects/single_cell_stitchit/work/scarlink/pbmc_input/top_1000_variable_genes.txt"

# shap genes
#scarlink_corr_file = "/projects/single_cell_stitchit/work/scarlink/pbmc_preprocessed_100pb_bins_test/model_performance/scarlink_test_spearman_corr_536.tsv"
#scarlink_corr_df = pd.read_csv(scarlink_corr_file, sep="\t", index_col=0)
#gene_ids = scarlink_corr_df.index
#performance_df = pd.read_csv(performance_file, sep="\t")
#performance_df = performance_df.fillna(0)
#top_2000_spearman = performance_df['Spearman'].nlargest(2000)
#gene_ids = top_2000_spearman.index.values
gene_ids = pd.read_csv(gene_file, sep="\t")["variable_genes"].values

#cell types from any feature matrix
fmat_path_0 = os.path.join(fmat_path, f"{gene_ids[0]}_lsi_activity_top_100.txt.gz")
sc_f_mat = pd.read_csv(fmat_path_0, sep="\t", index_col=0)
cell_ids = sc_f_mat.index

test_cells_file = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/test_cells.txt" #training cells as background
test_cells = np.loadtxt(test_cells_file, dtype=str)
test_idx = [cell_ids.tolist().index(x) for x in test_cells]
#test_idx = [item for item in test_idx] #if item not in bad_idx]

train_idx = [idx for idx in range(0, len(cell_ids)) if idx not in test_idx]
#train_idx = [item for item in train_idx] #if item not in bad_idx]

celltype_annot_path = "/projects/single_cell_stitchit/work/binned_RF/input/paired_PBMC10k/cell_type_annotation.tsv"
celltype_df = pd.read_csv(celltype_annot_path, sep="\t", index_col=0)
celltype_df =  celltype_df.reindex(cell_ids)
cell_types = celltype_df['celltypes'].unique()
celltype_df_train = celltype_df.iloc[train_idx, :] 

def explain_marginal(gene_model, fmat_train):
      explainer_tree = shap.TreeExplainer(gene_model, data =fmat_train, model_output="raw", feature_perturbation="interventional") #marginal distrubution: assumes independent features
      return explainer_tree

#treep-path dependent does not converge: trees too deep (depth ~500)
def explain_joint(gene_model):
      explainer_tree = shap.TreeExplainer(gene_model, data =None, model_output="raw", feature_perturbation="tree_path_dependent") #joint distribution: includes feature correlations, potential problem: high feature importance for feature that is not used by model, but has high correlation with feature that is important for the model 
      return explainer_tree

def compute_shap_values_global(explainer, x_train, cell_types, celltype_df_train, n_jobs):
    def global_shap_celltype(ct, x_train):
        ct_idx = (celltype_df_train['celltypes'] == ct).values
        #if np.sum(ct_idx) < 100:
            #return None  # Skip cell types with fewer than 100 samples
        x_ct_subset = x_train[ct_idx]
        shap_values_tree = explainer.shap_values(x_ct_subset)
        mean_shap = np.mean(shap_values_tree, axis=0)
        return ct, mean_shap
    shap_results = Parallel(n_jobs=n_jobs)(delayed(global_shap_celltype)(ct, x_train) for ct in np.unique(cell_types))
    shap_results = [res for res in shap_results if res is not None]
    ct_list, shap_list = zip(*shap_results)
    return ct_list, shap_list

#genes = []
#gene_id=gene_ids[0]
for gene_id in gene_ids:
    try:    
        print(gene_id)
        out_dir = f"{output_path}/{gene_id}/"
        if os.path.exists(out_dir):
            print(f"{out_dir} exists already. Skip gene.")
            continue
        model_file = os.path.join(model_path, f"models_{gene_id}_lsi_activity_top_100_sc_RF_0", f"{gene_id}_lsi_activity_top_100_sc_RF_0_model.joblib")
        print(model_file)
        if not os.path.exists(model_file):
            print(f"{model_file} does not exist. Skip gene.")
            continue
        #genes.append(gene_id)
        feature_mat_file = os.path.join(fmat_path, f"{gene_id}_lsi_activity_top_100.txt.gz")
        # expression prediction
        x = pd.read_csv(feature_mat_file, sep="\t", index_col=0)
        #TODO: log-transform feature matrix
        x = np.log2((1 + x.iloc[:, 0:(x.shape[1]-1)].astype(float)))
        x_train = x.iloc[train_idx, :] #log_transformed x_train
        #load gene model
        rf_model = joblib.load(model_file)
        # shap explainer model
        explainer = explain_marginal(rf_model,x_train)
        #shap_lst = [] 
        #zscore_d = {}
        #shap_d = {}
        ct_list, shap_lst = compute_shap_values_global(explainer, x_train, cell_types, celltype_df_train, n_jobs)
        shap_arr = np.vstack(shap_lst)
        z_scores = zscore(shap_arr, axis=1)
        #z_scores = zscore(np.ravel(shap_lst)).reshape(shap_lst.shape)    
        zscore_d = dict(zip(ct_list, z_scores))  
        shap_d = dict(zip(ct_list, shap_arr))

        for ct in cell_types:
            if ct not in zscore_d:
                #TODO: NA instead of 0?
                zscore_d[ct] = np.zeros(x.shape[1]) #for comparison with scarlink, shap values are not saved for cell-types with less than 100 cells to allow backtracking

        df_zscore = pd.DataFrame.from_dict(zscore_d, orient="index")
        df_zscore.columns = x.columns
        df_shap = pd.DataFrame.from_dict(shap_d, orient="index")  
        df_shap.columns = x.columns 

        os.makedirs(out_dir, exist_ok=True)

        df_zscore.to_csv(f"{out_dir}/{gene_id}_zscore_shap_sc_1MB500bp_q100.tsv", sep="\t", index=True, header=True)
        df_shap.to_csv(f"{out_dir}/{gene_id}_shap_sc_1MB500bp_q100.tsv", sep="\t", index=True, header=True)      

    except Exception as e:
        error_msg = f"Error processing {gene_id}: {str(e)}\n"
        print(error_msg)
        with open(log_file, "a") as f:
            f.write(error_msg)
