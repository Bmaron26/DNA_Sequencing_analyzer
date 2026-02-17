#!/bin/bash
#
# run_breseq.sh - Run breseq analysis on bacterial WGS samples
#
# Breseq is better than FreeBayes at detecting:
#   - Large deletions
#   - IS element insertions
#   - Gene duplications
#   - Complex structural variants
#
# Usage:
#   ./run_breseq.sh --reference <ref.gbk> --samples-dir <dir> --output-dir <dir>
#
# Or run individual samples:
#   ./run_breseq.sh --reference <ref.gbk> --r1 sample_1.fastq.gz --r2 sample_2.fastq.gz --name sample --output-dir <dir>
#
# Requirements:
#   - breseq (install via conda: conda install -c bioconda breseq)
#   - Reference in GenBank format (.gbk) for full annotation

set -e

# Default values
THREADS=4
MIN_COVERAGE=5

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --reference, -r    Reference genome (GenBank .gbk format recommended)"
    echo "  --samples-dir, -d  Directory containing FASTQ files (paired: *_1.fastq.gz, *_2.fastq.gz)"
    echo "  --output-dir, -o   Output directory for results"
    echo "  --r1               Forward reads (for single sample mode)"
    echo "  --r2               Reverse reads (for single sample mode)"
    echo "  --name, -n         Sample name (for single sample mode)"
    echo "  --threads, -t      Number of threads (default: 4)"
    echo "  --polymorphism     Enable polymorphism mode (for mixed populations)"
    echo "  --help, -h         Show this help message"
    echo ""
    echo "Examples:"
    echo "  # Batch mode - process all samples in a directory"
    echo "  $0 -r reference.gbk -d /path/to/fastq -o /path/to/output"
    echo ""
    echo "  # Single sample mode"
    echo "  $0 -r reference.gbk --r1 Mel1_1.fastq.gz --r2 Mel1_2.fastq.gz -n Mel1 -o results"
}

# Parse arguments
POSITIONAL_ARGS=()
POLYMORPHISM_MODE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --reference|-r)
            REFERENCE="$2"
            shift 2
            ;;
        --samples-dir|-d)
            SAMPLES_DIR="$2"
            shift 2
            ;;
        --output-dir|-o)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --r1)
            R1="$2"
            shift 2
            ;;
        --r2)
            R2="$2"
            shift 2
            ;;
        --name|-n)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --threads|-t)
            THREADS="$2"
            shift 2
            ;;
        --polymorphism)
            POLYMORPHISM_MODE="-p"
            shift
            ;;
        --help|-h)
            print_usage
            exit 0
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Check for breseq
if ! command -v breseq &> /dev/null; then
    echo -e "${RED}Error: breseq not found!${NC}"
    echo ""
    echo "Install breseq via conda:"
    echo "  conda install -c bioconda breseq"
    echo ""
    echo "Or via apt (Ubuntu/Debian):"
    echo "  sudo apt install breseq"
    exit 1
fi

# Validate inputs
if [ -z "$REFERENCE" ]; then
    echo -e "${RED}Error: Reference genome required (--reference)${NC}"
    print_usage
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo -e "${RED}Error: Reference file not found: $REFERENCE${NC}"
    exit 1
fi

if [ -z "$OUTPUT_DIR" ]; then
    echo -e "${RED}Error: Output directory required (--output-dir)${NC}"
    print_usage
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "======================================================================"
echo "                    BRESEQ ANALYSIS PIPELINE"
echo "======================================================================"
echo ""
echo "Reference: $REFERENCE"
echo "Output:    $OUTPUT_DIR"
echo "Threads:   $THREADS"
if [ -n "$POLYMORPHISM_MODE" ]; then
    echo "Mode:      Polymorphism (mixed population)"
else
    echo "Mode:      Consensus (clonal)"
fi
echo ""

# Function to run breseq on a single sample
run_breseq_sample() {
    local name=$1
    local r1=$2
    local r2=$3
    local outdir="$OUTPUT_DIR/$name"

    echo -e "${GREEN}Processing: $name${NC}"
    echo "  R1: $r1"
    echo "  R2: $r2"
    echo "  Output: $outdir"

    # Skip if already completed
    if [ -f "$outdir/output/index.html" ]; then
        echo -e "${YELLOW}  Skipping (already completed)${NC}"
        return 0
    fi

    # Run breseq
    breseq \
        -j $THREADS \
        -r "$REFERENCE" \
        $POLYMORPHISM_MODE \
        -o "$outdir" \
        "$r1" "$r2"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}  Completed successfully!${NC}"
    else
        echo -e "${RED}  Failed!${NC}"
        return 1
    fi
}

# Single sample mode
if [ -n "$R1" ] && [ -n "$R2" ] && [ -n "$SAMPLE_NAME" ]; then
    echo "Running single sample mode..."
    echo ""
    run_breseq_sample "$SAMPLE_NAME" "$R1" "$R2"

# Batch mode
elif [ -n "$SAMPLES_DIR" ]; then
    echo "Running batch mode..."
    echo "Scanning for samples in: $SAMPLES_DIR"
    echo ""

    # Find all paired FASTQ files
    declare -A SAMPLES

    for f in "$SAMPLES_DIR"/*_1.fastq.gz "$SAMPLES_DIR"/*_1.fq.gz "$SAMPLES_DIR"/*_R1.fastq.gz "$SAMPLES_DIR"/*_R1_001.fastq.gz; do
        if [ -f "$f" ]; then
            # Extract sample name
            basename=$(basename "$f")
            if [[ $basename =~ ^(.+)_1\.fastq\.gz$ ]]; then
                name="${BASH_REMATCH[1]}"
                r2="${f/_1.fastq.gz/_2.fastq.gz}"
            elif [[ $basename =~ ^(.+)_1\.fq\.gz$ ]]; then
                name="${BASH_REMATCH[1]}"
                r2="${f/_1.fq.gz/_2.fq.gz}"
            elif [[ $basename =~ ^(.+)_R1\.fastq\.gz$ ]]; then
                name="${BASH_REMATCH[1]}"
                r2="${f/_R1.fastq.gz/_R2.fastq.gz}"
            elif [[ $basename =~ ^(.+)_R1_001\.fastq\.gz$ ]]; then
                name="${BASH_REMATCH[1]}"
                r2="${f/_R1_001.fastq.gz/_R2_001.fastq.gz}"
            else
                continue
            fi

            if [ -f "$r2" ]; then
                SAMPLES[$name]="$f|$r2"
            fi
        fi
    done

    if [ ${#SAMPLES[@]} -eq 0 ]; then
        echo -e "${RED}No paired FASTQ files found!${NC}"
        echo "Expected naming: *_1.fastq.gz / *_2.fastq.gz"
        exit 1
    fi

    echo "Found ${#SAMPLES[@]} samples:"
    for name in "${!SAMPLES[@]}"; do
        echo "  - $name"
    done
    echo ""

    # Process each sample
    FAILED=()
    SUCCESS=()

    for name in "${!SAMPLES[@]}"; do
        IFS='|' read -r r1 r2 <<< "${SAMPLES[$name]}"

        if run_breseq_sample "$name" "$r1" "$r2"; then
            SUCCESS+=("$name")
        else
            FAILED+=("$name")
        fi
        echo ""
    done

    # Summary
    echo "======================================================================"
    echo "                         SUMMARY"
    echo "======================================================================"
    echo ""
    echo -e "${GREEN}Successful: ${#SUCCESS[@]}${NC}"
    for s in "${SUCCESS[@]}"; do
        echo "  - $s"
    done

    if [ ${#FAILED[@]} -gt 0 ]; then
        echo ""
        echo -e "${RED}Failed: ${#FAILED[@]}${NC}"
        for s in "${FAILED[@]}"; do
            echo "  - $s"
        done
    fi

else
    echo -e "${RED}Error: Either --samples-dir or (--r1, --r2, --name) required${NC}"
    print_usage
    exit 1
fi

echo ""
echo "======================================================================"
echo "                         RESULTS"
echo "======================================================================"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "For each sample, key output files:"
echo "  <sample>/output/index.html        - Interactive HTML report"
echo "  <sample>/output/output.gd         - GenomeDiff format (machine-readable)"
echo "  <sample>/output/evidence/         - Supporting evidence"
echo ""
echo "To view results, open the HTML file in a browser:"
echo "  firefox $OUTPUT_DIR/<sample>/output/index.html"
echo ""
echo "To compare with FreeBayes results, use:"
echo "  python scripts/compare_breseq_freebayes.py --breseq-dir $OUTPUT_DIR --freebayes-dir <freebayes_results>"
