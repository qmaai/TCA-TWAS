#! /bin/bash
# qsub -N lasso_after_tca -l h_data=24G,h_rt=03:00:00 -t 1-51:1 lasso_after_tca.sh -V
. /u/local/Modules/default/init/modules.sh
module load R/3.6.0

INDEX=${SGE_TASK_ID}
BULK=200
START=$(((INDEX-1)*BULK+1))
END=$((START+BULK-1))
stdfile="$START-$END.$(date "+%Y.%m.%d-%H.%M.%S").log"
rootfolder="/u/home/q/qmaai/project-sriram/tca-twas/"

Rscript $rootfolder/pipeline/lasso_after_tca.r --z_path $rootfolder/data/tca_result/full_sweep_without_too_many_cis_snps_gene/z_result_${START}_${END}.RData --mdl_path $rootfolder/data/tca_result/full_sweep_without_too_many_cis_snps_gene/mdl_result_${START}_${END}.RData --output_path $rootfolder/data/lasso_result/ > $rootfolder/data/lasso_result/log/$stdfile --verbose 2>&1

