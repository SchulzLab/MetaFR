import os
import signal
import threading
import gc
import numpy.ma as ma
import tracemalloc
import ctypes
import sys
import pickle
import gzip
import numpy as np
import math
import pandas as pd
import os
import sys
import statistics
import re
import joblib
import glob
###
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
###
from sklearn.linear_model import ElasticNet
from sklearn.ensemble import RandomForestRegressor
###
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
###
import scipy
from scipy import stats
#import scipy.stats

### Set a random seed for the results to be reproducible
import random

################################################################
################################################################
################################################################
################################################################
################################################################

def main():
    tracemalloc.start()

    ##### Read the data and begin constructing the model

    #TODO: paths models, results, plots
    #TODO: hyperparameters RandomForestRegressor() as user parameters in nextflow -> optimization

    #pre-processing: adapt input format: 1. read ATAC data from bigWig files (per metacell) and store as count matrix per gene: regions = equally sized bins in target gene window,
        # 2. add gene expression vector (response) from RNA count file, extract mean signal in bins? (lrumpf)
    #TODO: cross-validation procedure for performance evaluation? (lrumpf)


    args= sys.argv
    #args= ["", ] ## for testing and debugging
    print(args)
    fn= args[1] #input file with binned activity data
    model_type= args[2] #for now only RF is supported
    seed_value= int(args[3]) #random seed for reproducibility
    expr_threshold = float(args[4]) #threshold for variance of log2 transformed gene expression
    min_nonzero_expr = int(args[5]) #minimum percentage of non-zero expression counts
    min_mc_total = int(args[6]) #minimum number metacells required
    min_nonzero_feature_mat = int(args[7]) #minimum percentage non-zero feature values
    test_cells_path = args[8] #path to file with test metacell - nf pipeline train-test split in activity binning step
    plot_dir = args[9] #path to output directory for plots
    log_file_path = args[10] #path to log file
    models_out_dir = args[11] #path to output directory for model files
    performance_out_dir = args[12] #path to output directory for model performance files

    args = [arg if arg != 'None' else None for arg in sys.argv[1:]] 

    random.seed(seed_value)
    pattern = '[\w-]+?(?=\.)'
    fn_trim= re.search(pattern, fn).group()
    output_tag = fn_trim + "_" + str(seed_value)
    print(output_tag)

    ### create output directories if they do not exist ###
    if not os.path.exists(models_out_dir):
        os.makedirs(models_out_dir)
    if not os.path.exists(performance_out_dir):
        os.makedirs(performance_out_dir)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)    
    
    with PdfPages(f"{plot_dir}/results_{output_tag}.pdf") as pdf:
        try:

            file_name = fn_trim
            print(file_name)
            #file_name = "/projects/single_cell_stitchit/work/binned_RF/input/MPAL/mpal_binned_activity_ct_donor/ENSG00000134748_1MB100bp_binned.txt.gz"
            #file_name = "/projects/apog/work/IHEC/Binned_Activity/ENSG00000010610_BinnedActivity.txt.gz"
            gene_name= file_name.split("_")[0]

            data = pd.read_csv(fn, sep= "\t").to_numpy()
            row_names = data[:, 0]

            #TODO: variance filter genes + #non-zero entries, #non-zero entries metacells
            #lrumpf
            #cut-off sparsity gene expression vector: default min. 10% non-zero entries
            # filter constant expression vector
            
            #variance filter on log2 transformed gene expression
            log_expr = np.log2(1 + data[:, data.shape[1]-1].astype(float)) 
            expr_var = np.var(log_expr, dtype=np.float64)

            fig1, ax1 = plt.subplots(nrows=1, ncols=1)
            ax1.hist(log_expr, bins=15, color='skyblue', edgecolor='black')
            ax1.set_title(gene_name + ": expression distribution, variance=" + str(expr_var))
            ax1.set_xlabel("log2(expr+1)")
            ax1.set_ylabel("count")
            fig1.tight_layout()
            pdf.savefig(fig1)
            plt.close(fig1)


            if expr_threshold is not None:
            #filter all genes with variance 0    
                if expr_var == 0:
                    log_file = open(log_file_path,"a")
                    log_file.write(f"{gene_name}: Expression variance is 0 - failed gene\n")
                    log_file.close()
                    sys.exit(f"Expression variance is 0 for {gene_name}. Gene is skipped.")
                #filter genes with variance below threshold    
                if (expr_var < expr_threshold):
                    log_file = open(log_file_path,"a")
                    log_file.write(f"{gene_name}: Expression variance is less than {expr_threshold}\n")
                    log_file.close()
                sys.exit(f"Expression variance for {gene_name} is less than {expr_threshold}. Gene is skipped.")

            if min_nonzero_expr is not None:    
                if type(np.count_nonzero(data[:, data.shape[1]-1])) != int and len(np.count_nonzero(data[:, data.shape[1]-1])) < (min_nonzero_expr/100)*data[:, data.shape[1]-1].size:
                    log_file = open(log_file_path,"a")
                    log_file.write(f"{gene_name}: Expression vector has less than {min_nonzero_expr} non-zero entries - failed gene\n")
                    log_file.close()
                    sys.exit(f"Expression vector for {gene_name} has less than {min_nonzero_expr} non-zero entries. Gene is skipped.")

            #lrumpf idea: save genes by filtering metacells individually? problems downstream analysis? kick out genes with filtered metacells?
            #lrumpf
            #allow max. 90% zero-entries in feature matrix
            if min_nonzero_feature_mat is not None:
                if np.count_nonzero(data) < (min_nonzero_feature_mat/100)*(data.shape[0]*data.shape[1]):
                    log_file = open(log_file_path,"a")
                    log_file.write(f"{gene_name}: Number of non-zero entries in feature matrix is less than {min_nonzero_feature_mat}% - failed gene\n")
                    log_file.close()
                    sys.exit(f'Number of non-zero entries in feature matrix is less than {min_nonzero_feature_mat}%. Gene {gene_name} is skipped.')

            #lrumpf
            #skip gene completely if number metacells is < 30
            if min_mc_total is not None:
                if data.shape[0] < min_mc_total:
                    log_file = open(log_file_path,"a")
                    log_file.write(f"{gene_name}: Number of metacells is less than {min_mc_total} - failed gene\n")
                    log_file.close()
                    sys.exit('Number of metacells is less than {min_mc_total}. Gene {gene_name} is skipped.')

            # random sampling test samples
            # test_cells_num = math.floor((test_cells_percent / 100) * len(row_names))
            # test_cells = np.random.choice(row_names, size=test_cells_num, replace=False)
            test_cells = np.loadtxt(test_cells_path, dtype=str) # train/test split in activity binning step

            #train-test split before feature selection (in feature matrix generation)
            #test_cells = np.loadtxt(test_cells_path, dtype=str)
        
            y = np.log2(1 + data[:, data.shape[1]-1].astype(float))
            x = np.log2(1 + data[:, 1:(data.shape[1]-1)].astype(float))

            ## Split data into training and test sets
            test_idx = [row_names.tolist().index(x) for x in test_cells]
            test_idx= [item for item in test_idx] 

            train_idx = [idx for idx in range(0, len(y)) if idx not in test_idx]
            train_idx = [item for item in train_idx] 

            row_names_train= row_names[train_idx]
            row_names_test= row_names[test_idx]

            x_train = x[train_idx, :]
            y_train = y[train_idx]

            x_test = x[test_idx, :]
            y_test = y[test_idx]

            
        ########################################################################################################################################

            #TODO: per feature?
            #mn = min(y_train)
            #mx = max(y_train)

            #scaling training data
            #y_train = (y_train - mn) / (mx - mn);

            regr_rf = RandomForestRegressor(max_features= "sqrt", random_state= seed_value, oob_score= True, n_jobs = 1)
            regr_rf.fit(x_train, y_train)
            y_pred = regr_rf.predict(x_test)

            y_pred_train = regr_rf.predict(x_train)

            ## training error and cor
            train_cor = np.corrcoef(y_pred_train.squeeze(), y_train)[0, 1].round(3)
            #backscaling
            #y_pred_train = y_pred_train * (mx - mn) + mn
            #y_train = y_train * (mx - mn) + mn
            train_err = mean_squared_error(y_pred_train, y_train).round(8)
            ## test error and cor
            test_cor = np.corrcoef(y_pred.squeeze(),y_test)[0, 1].round(3)
            print(test_cor)
            test_spcor= stats.spearmanr(y_pred.squeeze(),y_test).statistic
            print(test_spcor)
            #test_cor = np.corrcoef(np.log2(1 + y_pred.squeeze()), np.log2(1 + y_test))[0, 1].round(2) # This is what Marcel suggested, which didn't make sense to me at all. I had to try it out nonetheless, and in fact this formulation made the results worse!

            #y_pred = y_pred * (mx - mn) + mn
            test_err = mean_squared_error(y_pred.squeeze(), y_test).round(8)

            print("RF test correlation")
            print(test_cor)

            print("variance of y_train= " + str(statistics.variance(y_train).round(5)))
            print("variance of y_test= " + str(statistics.variance(y_test).round(5)))
            
            fig2, ax2 = plt.subplots(nrows=1, ncols=1)
            ax2.scatter(y_pred, y_test)
            ax2.set_title(gene_name + " cor= " + str(test_cor) + ", error= "+ str(test_err))
            ax2.set_xlabel("RF Prediction")
            ax2.set_ylabel("Actual Expression")
            fig2.tight_layout()
            #plt.show()
            #plt.savefig('CNN_performance_' + gene_name + '.pdf')
            pdf.savefig(fig2)
            plt.close(fig2)

            #TODO: paths as user params nextflow (lrumpf)
            #if not os.path.exists('/projects/single_cell_stitchit/work/binned_RF/output/MPAL/models/models_' + output_tag):
            #    os.makedirs('/projects/single_cell_stitchit/work/binned_RF/output/MPAL/models/models_' + output_tag, exist_ok=True)
            #joblib.dump(regr_rf, '/projects/single_cell_stitchit/work/binned_RF/output/MPAL/models/models_' + output_tag  + '/' + output_tag + '_model.joblib')
            if not os.path.exists(models_out_dir + '/models_' + output_tag):
                os.makedirs(models_out_dir + '/models_' + output_tag, exist_ok=True)
            joblib.dump(regr_rf, models_out_dir + '/models_' + output_tag  + '/' + output_tag + '_model.joblib')
            print("######################################################################################################")

            #if not os.path.exists("/projects/single_cell_stitchit/work/binned_RF/output/MPAL/results/results_" + output_tag):
            #    os.makedirs("/projects/single_cell_stitchit/work/binned_RF/output/MPAL/results/results_" + output_tag, exist_ok=True)
            if not os.path.exists(performance_out_dir + "/results_" + output_tag):
                os.makedirs(performance_out_dir +  "/results_" + output_tag, exist_ok=True)

            df_cor = pd.DataFrame([train_cor, test_cor, test_spcor])
            df_cor.index= ["RF_train_cor", "RF_test_Pearson_cor", "RF_test_Spearman_cor"]
            df_err = pd.DataFrame([train_err, test_err])
            df_err.index= ["RF_train_err", "RF_test_err"]

            #df_cor.to_csv("/projects/single_cell_stitchit/work/binned_RF/output/MPAL/results/results_" + output_tag + "/cor_summary_w1MB_100bs.csv", sep= "\t")
            #df_err.to_csv("/projects/single_cell_stitchit/work/binned_RF/output/MPAL/results/results_" + output_tag + "/err_summary_w1MB_100bs.csv", sep= "\t")
            df_cor.to_csv(performance_out_dir + "/results_" + output_tag + "/cor_summary.csv", sep= "\t")
            df_err.to_csv(performance_out_dir + "/results/results_" + output_tag + "/err_summary.csv", sep= "\t")

        except Exception as e:
            with open(log_file_path, "a") as log_file:
                log_file.write(f"Error processing {gene_name}: {str(e)}\n")
            raise
        
        finally:
            libc = ctypes.CDLL("libc.so.6") #linux specific
            libc.malloc_trim(0)
            gc.collect()

if __name__ == "__main__":
    main()