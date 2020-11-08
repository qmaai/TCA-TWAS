#! /bin/bash
# qsub -N tca_gene_bulk -l h_data=24G,h_rt=03:00:00 -t 1-51:1 tca_gene_bulk.sh -V
. /u/local/Modules/default/init/modules.sh
module load R/3.6.0

INDEX=${SGE_TASK_ID}
BULK=200
START=$(((INDEX-1)*BULK+1))
END=$((START+BULK-1))
stdfile="$START-$END.$(date "+%Y.%m.%d-%H.%M.%S").log"
rootfolder="/u/home/q/qmaai/project-sriram/tca-twas/"
Rscript $rootfolder/pipeline/tca_gene_bulk.r --bulk_size $BULK --start_gene_pos $START > $rootfolder/data/tca_result/log/$stdfile 2>&1

