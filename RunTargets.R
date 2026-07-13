
# runs by using sb_run.sh to handle the coordination of workers -- everything in the shell file
# is all aobut the coordinator node.
# run the script by typing
# sbatch sb_run.sh
# in a terminal (not in R)

# set directories
setwd(this.path::this.dir())

# for cluster, cpp functions are compiled to these locations so it doesn't keep compiling them in parallel nodes
# if ('./cppcache_mv' %in% list.dirs() == FALSE){
#     dir.create('./cppcache_mv')
# }
# if ('./cppcache_ffoi' %in% list.dirs() == FALSE){
#     dir.create('./cppcache_ffoi')
# }
#
# # create space for log files in working directory
# if ('./log' %in% list.dirs() == FALSE){
#     dir.create('./log')
# }

#source targets file
source("_targets.R")

#Make pipeline
tar_make(callr_function = NULL, use_crew=FALSE, as_job=FALSE) # for troubleshooting
# tar_make_clustermq(workers = 8, log_worker = TRUE, reporter = 'verbose') # for parallelization/cluster implementation
# tar_make_clustermq(workers = 40, log_worker = TRUE, reporter = 'verbose') # for parallelization/cluster implementation

