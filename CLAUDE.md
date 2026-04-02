# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline project that combines Whole-Genome Sequencing (WGS) analysis with MitoEdit (MTE) for mitochondrial DNA variant analysis and TALEN design. The architecture consists of:

- **Frontend**: Streamlit UI (`UI Files/wgs_ui_v2.py`) for parameter collection
- **Backend**: FastAPI server (`api_server.py`) for job orchestration
- **Workflow Engine**: Nextflow DSL2 pipeline (`nf_integrated/nf_mte/mitoedit.nf`)
- **Analysis Tool**: MitoEdit for TALEN site prediction
- **Visualization**: IGV-Reports for interactive genome browsing

## Essential Commands

### Environment Setup
```bash
# Automated setup (installs Miniconda, Java 21, Nextflow v25.04.7, environments)
bash nf_integrated/getting-nf-conda-setup.sh

# Manual conda environment creation
cd nf_integrated/nf_mte
conda create -n mte -f mte.yml
conda create -n ui python=3.10
conda create -n wgs -y python=3.10 samtools=1.12 bwa=0.7.17
```

### Running the Pipeline Components

**Nextflow MitoEdit Pipeline:**
```bash
cd nf_integrated/nf_mte
conda activate mte
nextflow run mitoedit.nf \
    --mtdna_seq_path NC_012920.1.txt \
    --input_tab test_WGS_output.tab \
    --output_prefix test_output
```

**FastAPI Backend Server:**
```bash
conda activate ui
python api_server.py
# or
uvicorn api_server:app --host 0.0.0.0 --port 8000
```

**Streamlit UI (Local):**
```bash
cd "UI Files"
conda activate ui
streamlit run wgs_ui_v2.py
# Access at http://localhost:8501
```

**Streamlit UI (VM Deployment):**
```bash
screen -S streamlit
streamlit run wgs_ui_v2.py --server.port 8505
# Ctrl+A, D to detach
# Currently running at http://74.235.82.166:8505
```

## Architecture Overview

### Data Flow
```
User Input → Streamlit UI → FastAPI API → Nextflow Pipeline → MitoEdit Tool → Results
                ↓
            IGV-Reports Generation → Interactive Visualization
```

### Key Components

**API Server (`api_server.py`):**
- REST endpoints: `/api/pipeline/run`, `/api/pipeline/status/{job_id}`, `/api/pipeline/report/{job_id}`
- Background job execution using Nextflow
- In-memory job tracking (development only - needs Redis/database for production)

**Nextflow Pipeline (`mitoedit.nf`):**
- DSL2 workflow with parallel TALEN prediction
- Processes: `splitCsv` → `flatMap` → `MITOEDIT` → `COMBINE_RESULTS`
- Supports both local and SLURM execution profiles

**MitoEdit Integration:**
- CLI tool for TALEN binding site prediction
- Uses TAL effector binding code (HD→C, NI→A, NG→T, NH/NK→G)
- Parallel execution for multiple variant positions

### Input/Output Formats

**Input:**
- BAM files (tumor/germline pairs, indexed)
- Tab-separated WGS output with columns: `mtDNA_pos`, `Ref. Allele`, `Coverage`, `Percent_alternative_allele`
- Mitochondrial reference sequence (`NC_012920.1.txt`)

**Output:**
- `{output_prefix}_combined.tab` - merged results with TALEN predictions
- `wgs_igv_report.html` and `mito_igv_report.html` - interactive visualizations
- Individual TALEN output files: `talen_output_{position}_{base}.txt`

## Configuration Files

**Nextflow Config (`nf_integrated/nf_mte/nextflow.config`):**
- Pipeline parameters: `input_tab`, `mtdna_seq_path`, `min_spacer` (14), `max_spacer` (18)
- Execution profiles: `standard` (local), `slurm` (HPC)

**Conda Environment (`nf_integrated/nf_mte/mte.yml`):**
- Python 3.10.18, BioPython 1.85, Pandas 2.3.3, MitoEdit 1.0.2

## TAL Effector Design Rules

Critical domain knowledge for TALEN design (see `TAL_rules.md`):
- **Binding code**: HD→C, NI→A, NG→T, NN→G>A, NH→G (preferred), NK→G
- **Position 0**: Requires 5′ T for optimal binding
- **Mismatch tolerance**: 5′ end mismatches most detrimental, 3′ end more tolerated
- **Optimal repeat count**: ~15-20 repeats
- **TALEN spacer**: 12-21 bp between paired half-sites
- **Mitochondrial specificity**: Minimal chromatin, low CpG methylation

## Development Notes

### Current Status
- ✅ Functional Nextflow pipeline with parallel execution
- ✅ Streamlit UI for parameter collection
- ✅ FastAPI backend with job orchestration
- ✅ IGV-Reports integration
- ⚠️ In-memory job tracking (not production-ready)
- ❌ No automated test suite (manual testing with `test_WGS_output.tab`)

### Key Dependencies
- Java 21 JRE (required for Nextflow)
- Nextflow v25.04.7
- MitoEdit 1.0.2 (pip package)
- samtools 1.12, bwa 0.7.17
- ANNOVAR 2020-06-07, Bambino 1.0.jar

### Testing
Currently uses manual integration testing:
- Test data: `test_WGS_output.tab`, `NC_012920.1.txt`
- No formal test suite - run full pipeline with test data to validate

### Error Handling
- Nextflow processes use `errorStrategy = 'ignore'`
- API has basic error responses
- Production deployment needs comprehensive error handling and logging

## Project Structure

```
├── api_server.py              # FastAPI backend (300 lines)
├── UI Files/wgs_ui_v2.py     # Streamlit frontend (134 lines)
├── nf_integrated/nf_mte/     # Nextflow pipeline module
│   ├── mitoedit.nf           # Main workflow (146 lines)
│   ├── nextflow.config       # Pipeline configuration
│   ├── mte.yml              # Conda environment spec
│   └── test_WGS_output.tab  # Test input data
├── igv/                      # IGV-Reports configuration
├── TAL_rules.md             # Domain knowledge for TALEN design
└── README.md                # Project overview and task tracking
```