#!/bin/bash

INPUT_DIR="/data/WorkDir/Wei/AF3_local/AB56aC_EV1057"
OUTPUT_DIR="$INPUT_DIR/af_output"

MODEL_DIR="/data/AlphaFold/AF3_DB/model_parameter"
PUBLIC_DB_DIR="/data/AlphaFold/AF3_DB"

AF3_IMAGE="alphafold3:v3.0.4"

mkdir -p "$OUTPUT_DIR"

echo "Using AlphaFold 3 Docker image: $AF3_IMAGE"

for json_file in "$INPUT_DIR"/B56*.json; do

    echo "Processing $json_file..."

    docker run --rm \
        --volume "$INPUT_DIR":/root/af_input \
        --volume "$OUTPUT_DIR":/root/af_output \
        --volume "$MODEL_DIR":/root/models \
        --volume "$PUBLIC_DB_DIR":/root/public_databases \
        --gpus all \
        "$AF3_IMAGE" \
        python run_alphafold.py \
        --json_path="/root/af_input/$(basename "$json_file")" \
        --model_dir="/root/models" \
        --output_dir="/root/af_output/$(basename "$json_file" .json)" \
        --num_diffusion_samples=5

done
