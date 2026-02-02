# Windows Setup Guide for Bacterial Mutation Analyzer

This guide will help you set up and run the Bacterial Mutation Analyzer on Windows.

## Option 1: Windows Subsystem for Linux (WSL) - Recommended

WSL provides the best compatibility with bioinformatics tools.

### Step 1: Install WSL

1. Open **PowerShell as Administrator** (right-click PowerShell > Run as administrator)
2. Run this command:
   ```
   wsl --install
   ```
3. Restart your computer when prompted
4. After restart, Ubuntu will open and ask you to create a username and password

### Step 2: Install Required Tools in WSL

Open WSL (Ubuntu) and run these commands one by one:

```bash
# Update package lists
sudo apt update

# Install Python and pip
sudo apt install -y python3 python3-pip python3-venv

# Install bioinformatics tools
sudo apt install -y bwa samtools bcftools

# Optional: Install fastp for faster QC
sudo apt install -y fastp
```

### Step 3: Set Up the Analyzer

```bash
# Navigate to your project (adjust path as needed)
cd /mnt/c/Users/hayouka-lab/Documents/Bar/bacterial-mutation-analyzer

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install the package
pip install -e .
```

### Step 4: Run Analysis

Your Windows files are accessible under `/mnt/c/` in WSL. For example:
- `C:\Users\hayouka-lab\Documents` becomes `/mnt/c/Users/hayouka-lab/Documents`

```bash
# Activate environment (do this each time you open WSL)
source venv/bin/activate

# Run batch analysis
bma batch -e /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/samples.csv \
          -r /mnt/c/path/to/reference.fna \
          -a /mnt/c/path/to/annotation.gff3 \
          --species "S. aureus" \
          -o /mnt/c/Users/hayouka-lab/Documents/Bar/3rd_evolution/WGS/results
```

---

## Option 2: Native Windows (Limited)

Some tools have Windows versions, but this is less reliable.

### Step 1: Install Python

1. Download Python from https://www.python.org/downloads/
2. Run installer - **CHECK "Add Python to PATH"**
3. Click "Install Now"

### Step 2: Install the Analyzer

Open **Command Prompt** (cmd) and run:

```cmd
cd C:\Users\hayouka-lab\Documents\Bar\bacterial-mutation-analyzer
python -m venv venv
venv\Scripts\activate
pip install -e .
```

### Step 3: Install Bioinformatics Tools

You'll need to install Windows versions of:
- BWA: https://github.com/lh3/bwa (compile or find pre-built)
- SAMtools: http://www.htslib.org/download/
- BCFtools: http://www.htslib.org/download/

This is more complex - **WSL is recommended instead**.

---

## Quick Reference Commands

### Check if tools are installed:
```bash
bma check
```

### Run single sample analysis:
```bash
bma analyze -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
            -r reference.fna -a annotation.gff3 \
            -n sample_name -o results/
```

### Run batch analysis with multiple samples:
```bash
bma batch -e samples.csv -r reference.fna -a annotation.gff3 \
          --species "S. aureus" -o results/
```

### Run batch analysis with ancestral filtering:
```bash
# Include ancestral sample in your CSV, then use --ancestral flag
bma batch -e samples.csv -r reference.fna -a annotation.gff3 \
          --species "S. aureus" --ancestral "Anc1" -o results/
```

### Compare samples after analysis:
```bash
bma compare results/ -e samples.csv -o comparison/
```

---

## Your Files

Based on your setup, you should have:

**Reference files** (from your existing folder):
- `consensus.fna` or `consensus.fa` - Reference genome
- `consensus.gff3` - Gene annotations

**Sample files** (your 4 test samples):
- `Mel1_1.fastq.gz`, `Mel1_2.fastq.gz` - Melittin replicate 1
- `Mel2_1.fastq.gz`, `Mel2_2.fastq.gz` - Melittin replicate 2
- `Pex1_1.fastq.gz`, `Pex1_2.fastq.gz` - Pexiganan replicate 1
- `Pex2_1.fastq.gz`, `Pex2_2.fastq.gz` - Pexiganan replicate 2

---

## Troubleshooting

### "Command not found" errors
- Make sure you activated the virtual environment: `source venv/bin/activate`

### "Tool not available" errors
- Run `bma check` to see which tools are missing
- Install missing tools with `sudo apt install <tool-name>`

### File path issues in WSL
- Windows paths need to be converted: `C:\folder` becomes `/mnt/c/folder`
- Use forward slashes `/` not backslashes `\`

### Permission denied
- Run with `sudo` if needed
- Or change file permissions: `chmod +x filename`
