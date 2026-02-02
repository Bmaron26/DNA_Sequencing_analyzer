#!/bin/bash
#
# Bacterial Mutation Analyzer - Quick Run Script
#
# Usage: ./run_analysis.sh <samples.csv> <reference.fna> <annotation.gff3> <output_dir>
#
# Example:
#   ./run_analysis.sh samples.csv reference.fna genes.gff3 results/
#

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Bacterial Mutation Analyzer${NC}"
echo -e "${CYAN}========================================${NC}"

# Check arguments
if [ $# -lt 4 ]; then
    echo -e "${RED}Usage: $0 <samples.csv> <reference.fna> <annotation.gff3> <output_dir> [ancestral_sample_id]${NC}"
    echo ""
    echo "Arguments:"
    echo "  samples.csv        - CSV file with sample information"
    echo "  reference.fna      - Reference genome FASTA file"
    echo "  annotation.gff3    - Gene annotation file (GFF3 format)"
    echo "  output_dir         - Directory for results"
    echo "  ancestral_sample   - (Optional) Sample ID of ancestral/WT strain"
    echo ""
    echo "Example:"
    echo "  $0 samples.csv consensus.fna consensus.gff3 results/"
    echo "  $0 samples.csv consensus.fna consensus.gff3 results/ Anc1"
    exit 1
fi

SAMPLES=$1
REFERENCE=$2
ANNOTATION=$3
OUTPUT=$4
ANCESTRAL=${5:-""}

# Check if files exist
if [ ! -f "$SAMPLES" ]; then
    echo -e "${RED}Error: Sample file not found: $SAMPLES${NC}"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo -e "${RED}Error: Reference file not found: $REFERENCE${NC}"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo -e "${RED}Error: Annotation file not found: $ANNOTATION${NC}"
    exit 1
fi

# Check if bma command is available
if ! command -v bma &> /dev/null; then
    echo -e "${RED}Error: 'bma' command not found${NC}"
    echo "Please activate the virtual environment first:"
    echo "  source venv/bin/activate"
    exit 1
fi

# Check tools
echo -e "${CYAN}Checking required tools...${NC}"
bma check

if [ $? -ne 0 ]; then
    echo -e "${RED}Error: Some required tools are missing${NC}"
    echo "Please install them before running the analysis"
    exit 1
fi

echo ""
echo -e "${GREEN}Starting analysis...${NC}"
echo "  Samples:    $SAMPLES"
echo "  Reference:  $REFERENCE"
echo "  Annotation: $ANNOTATION"
echo "  Output:     $OUTPUT"

if [ -n "$ANCESTRAL" ]; then
    echo "  Ancestral:  $ANCESTRAL (will filter pre-existing variants)"
fi

echo ""

# Build command
CMD="bma batch -e $SAMPLES -r $REFERENCE -a $ANNOTATION -o $OUTPUT --species 'S. aureus' -v"

if [ -n "$ANCESTRAL" ]; then
    CMD="$CMD --ancestral $ANCESTRAL"
fi

# Run analysis
echo -e "${CYAN}Running: $CMD${NC}"
echo ""

eval $CMD

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Analysis Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo "Results saved to: $OUTPUT"
    echo ""
    echo "Key output files:"
    echo "  - comparison_report.json    : Multi-sample comparison"
    echo "  - mutation_matrix.csv       : Sample x mutation matrix"
    echo "  - convergent_mutations.csv  : Mutations in multiple samples"

    if [ -n "$ANCESTRAL" ]; then
        echo "  - novel_mutations_report.json : Mutations NOT in ancestral"
        echo "  - novel_mutation_matrix.csv   : Novel mutations only"
    fi

    echo ""
    echo "Each sample folder contains:"
    echo "  - {sample}_mutations.csv    : All mutations"
    echo "  - {sample}_amr_report.json  : AMR database annotations"
else
    echo ""
    echo -e "${RED}Analysis failed. Check the error messages above.${NC}"
    exit 1
fi
