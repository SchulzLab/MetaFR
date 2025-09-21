import scarlink as sl
import scarlink.src.visualization as scv
import pandas
import subprocess

#TODO: run on GPU?

#outdir = "pbmc_test"
#outdir = "/projects/single_cell_stitchit/work/scarlink/pbmc_output_100pb_bins"
outdir="/projects/single_cell_stitchit/work/scarlink/pbmc_output_500pb_bins_train"
#subprocess.run(["scarlink_processing --scrna /projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_correct_cbs.rds --scatac /projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scatac -o /projects/single_cell_stitchit/work/scarlink/pbmc_output_100pb_bins -nc 25"], shell=True, capture_output=True, text=True)
#subprocess.run(["time scarlink -o /projects/single_cell_stitchit/work/scarlink/pbmc_output_100pb_bins -g hg38 -np 25"], shell=True, capture_output=True, text=True)

#initial run: run without SHAP calculation
subprocess.run(["time scarlink_processing --scrna /projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_training_cells_celltype_annotation_cell_prefix_500bp.rds --scatac /projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_train_500bp -o /projects/single_cell_stitchit/work/scarlink/pbmc_output_500pb_bins_train -nc 25"], shell=True, capture_output=True, text=True)
subprocess.run(["time scarlink -o /projects/single_cell_stitchit/work/scarlink/pbmc_output_500pb_bins_train -g hg38 -np 25"], shell=True, capture_output=True, text=True)

######################################################################################################################################################################################################################################

# pre-processing test cells
subprocess.run(["time scarlink_processing --scrna /projects/single_cell_stitchit/work/scarlink/pbmc_input/pbmc_scrna_test_cells_cell_prefix_500bp.rds --scatac /projects/single_cell_stitchit/work/scarlink/pbmc_input/250kb_500bp_setup/pbmc_scatac_test_500bp -o /projects/single_cell_stitchit/work/scarlink/pbmc_preprocessed_500pb_bins_test -nc 25"], shell=True, capture_output=True, text=True)


