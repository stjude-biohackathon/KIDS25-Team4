# WGS Pipeline UI

Streamlit interface for WGS mtDNA variant analysis parameter collection.

## wgs_ui_v2.py
- Reference genome: hg38 only
- Pipeline selection: WGS only / WGS + MitoEdit (full pipeline)
- SNV class filter: TU (Tumor-enriched), HS (Shared), GH (Germline heteroplasmy), INH (Inherited)
- VAF threshold: Restricted to >0 and <99.9
- Access at: http://74.235.82.166:8505

## Installation & Setup

### Requirements
- Python 3.10+
- Streamlit

### Local Setup Instructions

1. Create a conda environment:  
`conda create -n ui python=3.10`  
`conda activate ui`  
`pip install streamlit`  

2. Run the application:  
`streamlit run wgs_ui_v2.py`

The interface will be accessible at http://localhost:8501 by default.

### Deployment on VM

The UI is currently running persistently on the team VM using a terminal multiplexer (screen).

**To start in screen:**  
`screen -S streamlit`  
`streamlit run wgs_ui_v2.py --server.maxUploadSize 51240 --server.port 8505`

Press Ctrl+A, then D to detach

**To reconnect to the session:**

`screen -r streamlit`  

**To list all screen sessions:**
`screen -ls`  

This ensures the UI remains accessible even when terminal sessions are closed or SSH connections are lost.

## Parameters Collected
- Paired BAM files (tumor/germline)
- Project name
- Output directory
- Reference genome
- Pipeline selection
- SNV class filter
- VAF threshold
- Minimum coverage
- Advanced filtering options

## Backend Integration

