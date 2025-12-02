import streamlit as st
import json
import pandas as pd

# Page configuration
st.set_page_config(page_title="WGS mtDNA Pipeline", layout="wide")

# Title
st.title("WGS mtDNA Variant Analysis Pipeline")
st.markdown("Upload your BAM files and configure analysis parameters")

# Sidebar for parameters
st.sidebar.header("Analysis Parameters")

# Main content - File uploads
st.header("1. Upload BAM Files")
col1, col2 = st.columns(2)

with col1:
    tumor_bam = st.file_uploader("Tumor BAM File", type=['bam'], help="Upload the tumor sample BAM file")
    
with col2:
    germline_bam = st.file_uploader("Germline/Normal BAM File", type=['bam'], help="Upload the matched normal sample BAM file")

# Project Information
st.header("2. Project Information")
col3, col4 = st.columns(2)

with col3:
    project_name = st.text_input("Project Name", placeholder="e.g., Patient_001", help="Unique identifier for this analysis")
    
with col4:
    output_dir = st.text_input("Output Directory", value="results/", help="Path where results will be saved")

# Reference Genome - hg38 only
st.header("3. Reference Genome")
reference_genome = st.selectbox(
    "Select Reference Genome",
    ["hg38"],
    help="Human genome reference version (only hg38 supported)"
)

# Pipeline Selection
st.header("4. Pipeline Selection")
pipeline_choice = st.selectbox(
    "Select Pipeline to Run",
    ["WGS only", "WGS + MitoEdit (full pipeline)"],
    help="Choose whether to run only variant calling or the full pipeline with base editing analysis"
)

# SNV Class Filter
st.header("5. SNV Class Filter")
snv_class_filter = st.selectbox(
    "Filter for SNV Class",
    ["TU (Tumor-enriched)", "HS (Shared)", "GH (Germline heteroplasmy)", "INH (Inherited)"],
    help="Select which SNV class to analyze"
)

# Variant Filters
st.header("6. Variant Filtering Parameters")
col5, col6 = st.columns(2)

with col5:
    min_vaf = st.number_input(
        "Minimum VAF Threshold (%)",
        min_value=0.1,
        max_value=99.9,
        value=3.0,
        step=0.1,
        help="Minimum variant allele fraction to report (must be >0 and <99.9)"
    )

with col6:
    min_coverage = st.number_input(
        "Minimum Coverage",
        min_value=10,
        max_value=10000,
        value=100,
        help="Minimum read depth required"
    )

# Advanced options (collapsible)
with st.expander("Advanced Options"):
    strand_bias_filter = st.checkbox("Filter strand bias", value=True)
    remove_common_variants = st.checkbox("Remove common variants", value=False)

# Submit button
st.markdown("---")
if st.button("Submit Analysis", type="primary"):
    # Validation
    errors = []
    
    if not tumor_bam:
        errors.append("Tumor BAM file is required")
    if not germline_bam:
        errors.append("Germline BAM file is required")
    if not project_name:
        errors.append("Project name is required")
    if not output_dir:
        errors.append("Output directory is required")
    if min_vaf <= 0 or min_vaf >= 99.9:
        errors.append("VAF threshold must be >0 and <99.9")
    
    if errors:
        st.error("Please fix the following errors:")
        for error in errors:
            st.write(f"- {error}")
    else:
        # Collect all parameters
        parameters = {
            "tumor_bam": tumor_bam.name if tumor_bam else None,
            "germline_bam": germline_bam.name if germline_bam else None,
            "project_name": project_name,
            "output_dir": output_dir,
            "reference_genome": reference_genome,
            "pipeline_choice": pipeline_choice,
            "snv_class_filter": snv_class_filter,
            "min_vaf": min_vaf / 100,  # Convert percentage to decimal
            "min_coverage": min_coverage,
            "strand_bias_filter": strand_bias_filter,
            "remove_common_variants": remove_common_variants
        }
        
        st.success("All inputs validated successfully!")
        
        # Display collected parameters
        st.subheader("Collected Parameters")
        st.json(parameters)
        
        # This is where backend integration would happen
        st.info("These parameters are ready to be passed to the Nextflow pipeline")
        
        # Write all parameters to a JSON file for backend processing
        with open("wgs_mito_pipeline_parameters.json", "w") as json_file:
            json.dump(parameters, json_file, indent=2)
        
        # Create a CSV sample sheet named Mit_WGS_SampleSheet.csv
        sample_sheet = pd.DataFrame({
            "Sample_Name": [project_name, project_name],
            "Sample_Type": ["Normal", "Tumor"],
            "BAM": [germline_bam.name if germline_bam else "", \
                tumor_bam.name if tumor_bam else ""]
        })
        sample_sheet.to_csv("Mit_WGS_SampleSheet.csv", index=False, sep="\t", header=False)
