# update conda
conda update -n base -c defaults conda

# create environment for running the snakemake pipeline
conda create --name snakemake -c bioconda -y snakemake

# setup gocommands
GOCMD_VER=$(curl -L -s https://raw.githubusercontent.com/cyverse/gocommands/main/VERSION.txt); \
curl -L -s https://github.com/cyverse/gocommands/releases/download/${GOCMD_VER}/gocmd-${GOCMD_VER}-darwin-amd64.tar.gz | tar zxvf -

# copy signalP to working directory
gocmd get --progress /iplant/home/shared/NCEMS/working-groups/energetic-origins/signalp6_fast 2-sequence-parsing/python-scripts/

# copy the signaP model files to the correct directory
cp 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/models/distilled_model_signalp6.pt 2-sequence-parsing/python-scripts/signalp6_fast/signalp-6-package/signalp/model_weights/

# activate the snakemake pipeline
conda activate snakemake
