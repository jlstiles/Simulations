########################################
# General setup

# Directory where sbatch-r.sh, sbatch-rmd.sh, etc. can be found.
#SCRIPT_DIR=scripts
SCRIPT_DIR=.

# Directory to store command results; set to "." to be current directory.
# OUTPUT_DIR=output
OUTPUT_DIR=.

# How do we want to run tasks? Can be slurm or bash currently.
# Use SLURM if possible, otherwise use bash.
# Can override if desired: "export JOB_ENGINE=shell"
ifndef JOB_ENGINE
  # Detect if we can use slurm, otherwise use shell.
  ifeq (, $(shell which sbatch))
		JOB_ENGINE=shell
	else
		JOB_ENGINE=slurm
	endif
	# TODO: check for SGE.
endif

######
# Savio configuration.

# This allows us to use environmental variables to override this default.
# e.g. we run in BASH: "export ACCOUNT=co_otheraccount"
ifndef ACCOUNT
	ACCOUNT=co_biostat
endif

# This allows us to use environmental variables to override this default.
ifndef PARTITION
	PARTITION=savio2
endif

# This allows us to override the default QOS by setting an environmental variable.
# e.g. we run in BASH: "export QOS=biostat_normal"
ifndef QOS
	# Choose one QOS and comment out the other, or use environmental variables.
	QOS=biostat_savio2_normal
	#QOS=savio_lowprio
endif

########################################
# Execution engines.

# Sbatch runs a SLURM job, e.g. on Savio or XSEDE.
SBATCH=sbatch -A ${ACCOUNT} -p ${PARTITION} --qos ${QOS}

# Setup R to run commands in the background and keep running after logout.
R=nohup nice -n 19 R CMD BATCH --no-restore --no-save

# TODO: support Sun Grid Engine (SGE) for grizzlybear2.
# Or just convert to batchtools?

########################################
# Tasks that can be run.

# Run an R file via "make analysis"
analysisTest: test.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

analysisHalvsDelta: caseHal.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

analysisbbd: master_script1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif


simkaraYMmis: sim_kara_YMmis.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYZmis: sim_kara_YZmis.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYSmis: sim_kara_YSmis.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkarawell: sim_kara_well.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYAwell: sim_kara_YAwell.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYmis: sim_kara_Ymis.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaracrazy: sim_kara_crazy.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYmistest: sim_kara_Ymis_test.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkarawelltest: sim_kara_well_test.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

clustertest: cluster_test.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYMmis1: sim_kara_YMmis1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYZmis1: sim_kara_YZmis1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYSmis1: sim_kara_YSmis1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkarawell1: sim_kara_well1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYAwell1: sim_kara_YAwell1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYmis1: sim_kara_Ymis1.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif


simkaraYMmis2: sim_kara_YMmis2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYZmis2: sim_kara_YZmis2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYSmis2: sim_kara_YSmis2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkarawell2: sim_kara_well2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYAwell2: sim_kara_YAwell2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

simkaraYmis2: sim_kara_Ymis2.R
ifeq (${JOB_ENGINE},slurm)
	${SBATCH} --nodes 1 --cpus-per-task=24 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif



# Options customized based on "7. GPU job script" at:
# http://research-it.berkeley.edu/services/high-performance-computing/running-your-jobs
gpu-test: gpu-test.Rmd
	sbatch -A ${ACCOUNT} -p savio2_gpu --qos savio_lowprio --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-rmd.sh --file=$< --dir=${OUTPUT_DIR}

# Launch a bash session on 2 compute nodes for up to 12 hours via "make bash".
bash:
	srun -A ${ACCOUNT} -p ${PARTITION}  -N 2 -t 12:00:00 --pty bash

####
# Add other rules here.
####

# Clean up any logs or temporary files via "make clean"
# Next line ensures that this rule works even if there's a file named "clean".
.PHONY : clean
clean:
	rm -f *.Rout
	rm -f slurm*.out
	rm -f install*.out
	rm -f cache/*
