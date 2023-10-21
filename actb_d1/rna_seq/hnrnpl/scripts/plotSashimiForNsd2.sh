#!/bin/bash
#SBATCH --account=rpp-jfcote11
#SBATCH --mail-user=poorani.subramani@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000
#SBATCH --time=3:00:0
#SBATCH --tmp=100G
#SBATCH --output="./slurm_outputs/%x_slurm-%j.out"

set -e
set -u
set -o pipefail

SLRM_LOCAL_TMP_DIR=$SLURM_TMPDIR/$RANDOM/
mkdir -p ${SLRM_LOCAL_TMP_DIR}

#---------------------------------------------------------#
module load gcc/9.3.0
module load gsl/2.6
module load cmake/3.23.1
module load python/3.8.10

# getting Nsd2 rMATS events if present
echo $(date -Iseconds)" SPG_RNAseq_pipeline: 8. Getting Nsd2 rMATS events if any.."
Rscript ${SCRIPT_CUR_DIR}/getNsd2rMATSevent.R
source ${DATA_INFO_DIR}/pipeline_config_vars.sh

#---------------------------------------------------------#

# if Sirt rMATS event is present, plot sashimi plots with rMATS2sashimiplot
if [[ $NSD2_RMATS_EVENT == "Yes" ]]
then

module load python/3.8.10
virtualenv --no-download $SLURM_TMPDIR/env
set +u
source $SLURM_TMPDIR/env/bin/activate
set -u
pip install --no-index --upgrade pip
pip install --no-index -r ${SCRIPT_CUR_DIR}/rMATS2sashimi_requirements.txt
module load samtools/1.16.1
module load bamtools/2.5.1

RMATS2SASHIMI_CODE_DIR="/home/subrampg/local_installs/rmats2sashimiplot/src/rmats2sashimiplot/"

python $RMATS2SASHIMI_CODE_DIR/rmats2sashimiplot.py \
	-o ${RMATS_DIR}/rMATS2sashimiPlot_outputs/ \
	--l1 $REFR_GROUP_NAME --l2 $TEST_GROUP_NAME \
	--event-type SE -e ${RMATS_DIR}/Nsd2_rMATS_SE_event.txt \
	--b1 ${RMATS_DIR}/input_refr_bamFiles.txt \
	--b2 ${RMATS_DIR}/input_test_bamFiles.txt \
	--group-info ${RMATS_DIR}/groupInfoForSashimi.gf \
	--fig-width 1.5 --fig-height 4 \
	--font-size 5 --color '#807D7D,#CC0011' \
	--intron_s 15 --no-text-background

if [[ $GENOME == "mm10" ]]
then
	COORDS_VAL="chr10:-:63320000:63339250:$GFF3_FILE"
elif [[ $GENOME == "hg38" ]]
then
	COORDS_VAL="chr10:+:67884650:67918000:$GFF3_FILE"
fi

python $RMATS2SASHIMI_CODE_DIR/rmats2sashimiplot.py \
	-o ${RMATS_DIR}/rMATS2sashimiPlot_outputs/ \
	--l1 $REFR_GROUP_NAME --l2 $TEST_GROUP_NAME \
	-c ${COORDS_VAL} \
	--b1 ${RMATS_DIR}/input_refr_bamFiles.txt \
	--b2 ${RMATS_DIR}/input_test_bamFiles.txt \
	--group-info ${RMATS_DIR}/groupInfoForSashimi.gf \
	--fig-width 3.5 --fig-height 4 \
	--font-size 5 --color '#807D7D,#CC0011' \
	--intron_s 15 --no-text-background

deactivate
rm -rf $SLURM_TMPDIR/env

fi

#---------------------------------------------------------#
#making symlinks
for SYML_DIR in RMATS_SYML_DIR
do
        SCRATCH_DIR_VARNAME="$(sed "s:_SYML::" <(echo ${SYML_DIR}))"
        eval SCRATCH_DIR='$'"${SCRATCH_DIR_VARNAME}"
        echo $SCRATCH_DIR_VARNAME" = "$SCRATCH_DIR
        for file in $(ls ${SCRATCH_DIR})
        do
        if [ ! -L "${!SYML_DIR}${file}" ]; then
                ln -s ${SCRATCH_DIR}${file} ${!SYML_DIR}${file}
        fi
        done

find ${!SYML_DIR} -xtype l -delete

done

