step 1: login to osc HPCs (cardinal)

ssh {your_user_name}@cardinal.osc.edu

step 2: prepare json file

{
    "name": "E06_POVPC",
    "modelSeeds": [1,2,3,4,5],
    "sequences": [
        {"protein": {
            "id": ["A"],
            "sequence": "MTQSPTFLAVTASKKVTISCTASESLYSSKHKVHYLAWYQKKPEQSPKLLIYGASNRYIGVPDRFTGSGSGTDFTLTISSVQVEDLTHYYCAQFYSYPLTFGAGTKLELK"
            }
        },
       {"protein": {
            "id": ["B"],
            "sequence": "MKLWLNWVFLLTLLHGIQCEVKLVESGGGLVQPGGSLRLSCATSGFTFSDFYMEWVRQPPGKRLEWIAASRNKANDYTTEYSASVKGRFIVSRDTSQSILYLQMNALRAEDTAIYYCARDYYGSSYWYFDVWGAGTTVTVSSES"
            }
        },
        {"ligand": {
            "id": ["C"],
            "smiles": "CCCCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCC=O"
            }
        }
    ],
    "dialect": "alphafold3",
    "version": 3
}

step 3: prepare slurm file

#!/bin/bash
#SBATCH --job-name=bNbs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --gpus-per-node=1
#SBATCH --account=PDS0370
#SBATCH --time=12:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wxh180@case.edu
#SBATCH --mem=128G
#SBATCH --output=logs/bNbs_%j.out
#SBATCH --error=logs/bNbs_%j.err

# Ensure the script runs in Bash
if [ -z "$BASH_VERSION" ]; then
    echo "Error: This script must be run with Bash." >&2
    exit 1
fi

# Reset module environment
module reset
module load alphafold3/3.0.1

# Define variables
MODEL_SRC="/users/PDS0295/cwr0487/AF3/parameters"
JSON_DIR="$(pwd -P)"
WORKDIR="/fs/scratch/PDS0370/AF3"
OUTPUT_DIR="${WORKDIR}/output/lipid_mAb"

# Debug: Print variable values
echo "Debug: MODEL_SRC=${MODEL_SRC}" >&2
echo "Debug: JSON_DIR=${JSON_DIR}" >&2
echo "Debug: WORKDIR=${WORKDIR}" >&2
echo "Debug: OUTPUT_DIR=${OUTPUT_DIR}" >&2

# Verify paths exist and are accessible
if [ ! -d "${MODEL_SRC}" ]; then
    echo "Error: Model directory ${MODEL_SRC} does not exist or is not accessible." >&2
    exit 1
fi

# Create and verify output directory
mkdir -p "${OUTPUT_DIR}"
if [ ! -w "${OUTPUT_DIR}" ]; then
    echo "Error: Output directory ${OUTPUT_DIR} is not writable." >&2
    exit 1
fi

# Verify run_alphafold.sh is available
if ! command -v run_alphafold.sh &> /dev/null; then
    echo "Error: run_alphafold.sh not found in PATH after loading alphafold3/3.0.1 module." >&2
    exit 1
fi

# Run AlphaFold 3
#echo "Running AlphaFold 3..." >&2
#run_alphafold.sh --model_dir="${MODEL_SRC}" \
#                 --output_dir="${OUTPUT_DIR}" \
#                 --json_path="${JSON_DIR}/${fn}.json"

# Iterate through all .json files in JSON_DIR
for json_file in "${JSON_DIR}"/*.json; do
    # Check if there are any JSON files
    if [[ -f "${json_file}" ]]; then
        # Extract the base filename without the .json extension
        fn=$(basename "${json_file}" .json)

        # Print status message to stderr
        echo "Running AlphaFold 3 for ${fn}..." >&2

        # Execute the AlphaFold command
        run_alphafold.sh --model_dir="${MODEL_SRC}" \
                         --output_dir="${OUTPUT_DIR}" \
                         --json_path="${json_file}"
    else
        echo "No JSON files found in ${JSON_DIR}" >&2
        exit 1
    fi
done

if [ $? -ne 0 ]; then
    echo "Error: run_alphafold.sh failed." >&2
    exit 1
fi
