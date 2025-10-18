"""
Streamlit GUI for Codon Usage Analysis
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import io

# Handle both relative and absolute imports
try:
    from .core.analysis import (
        parse_fasta_from_text,
        compute_codon_frequencies,
        compute_rscu_weights,
        compute_amino_acid_usage,
        compute_per_gene_rscu,
        compute_codon_usage_per_1000,
        compute_relative_codon_frequencies
    )
except ImportError:
    from codon_usage_gui.core.analysis import (
        parse_fasta_from_text,
        compute_codon_frequencies,
        compute_rscu_weights,
        compute_amino_acid_usage,
        compute_per_gene_rscu,
        compute_codon_usage_per_1000,
        compute_relative_codon_frequencies
    )


def create_rscu_heatmap(df_rscu):
    """Create a heatmap of RSCU values."""
    # Pivot the data for heatmap
    pivot_df = df_rscu.pivot_table(
        index='Amino_Acid',
        columns='Codon',
        values='RSCU',
        fill_value=0
    )

    fig = px.imshow(
        pivot_df,
        aspect="auto",
        color_continuous_scale="RdYlBu_r",
        title="Relative Synonymous Codon Usage (RSCU) Heatmap"
    )
    fig.update_layout(
        xaxis_title="Codons",
        yaxis_title="Amino Acids",
        height=600
    )
    return fig


def create_codon_usage_bar_plot(df_codcount):
    """Create a bar plot of codon usage frequencies."""
    # Remove STOP codons for cleaner visualization
    df_plot = df_codcount[df_codcount['Amino_Acid'] != 'STOP'].copy()

    fig = px.bar(
        df_plot,
        x='Codon',
        y='Obs_Freq',
        color='Amino_Acid',
        title="Codon Usage Frequencies",
        labels={'Obs_Freq': 'Frequency', 'Codon': 'Codon'}
    )
    fig.update_layout(
        xaxis_title="Codons",
        yaxis_title="Frequency",
        height=500,
        xaxis={'categoryorder': 'total descending'}
    )
    return fig


def create_aa_usage_comparison(aa_df):
    """Create a comparison plot of expected vs observed amino acid usage."""
    fig = go.Figure()

    fig.add_trace(go.Bar(
        name='Expected',
        x=aa_df['Amino_acid'],
        y=aa_df['Expected_Freq(%)'],
        marker_color='lightblue'
    ))

    fig.add_trace(go.Bar(
        name='Observed',
        x=aa_df['Amino_acid'],
        y=aa_df['Obs_Freq(%)'],
        marker_color='darkred'
    ))

    fig.update_layout(
        title="Expected vs Observed Amino Acid Usage",
        xaxis_title="Amino Acids",
        yaxis_title="Frequency (%)",
        barmode='group',
        height=500
    )
    return fig


def create_rscu_distribution(df_rscu):
    """Create a histogram of RSCU distribution."""
    fig = px.histogram(
        df_rscu[df_rscu['Amino_Acid'] != 'STOP'],
        x='RSCU',
        nbins=30,
        title="Distribution of RSCU Values",
        labels={'RSCU': 'RSCU Value', 'count': 'Number of Codons'}
    )
    fig.add_vline(x=1, line_dash="dash", line_color="red",
                  annotation_text="Expected value (1.0)")
    return fig


def main():
    st.set_page_config(
        page_title="Codon Usage Analysis",
        page_icon="🧬",
        layout="wide"
    )

    st.title("🧬 Codon Usage Bias Analysis Tool")
    st.markdown("""
    This tool provides comprehensive analysis of codon usage bias from FASTA sequences.
    Upload your coding sequences to analyze RSCU, amino acid usage, and codon frequencies.
    """)

    # Sidebar for input
    st.sidebar.header("📁 Input Data")

    input_method = st.sidebar.radio(
        "Choose input method:",
        ["Upload FASTA file", "Paste FASTA text"]
    )

    sequences = None
    headers = None

    if input_method == "Upload FASTA file":
        uploaded_file = st.sidebar.file_uploader(
            "Choose a FASTA file",
            type=['fasta', 'fa', 'fas', 'txt'],
            help="Upload a FASTA file containing coding sequences"
        )

        if uploaded_file is not None:
            stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            fasta_text = stringio.read()
            headers, sequences = parse_fasta_from_text(fasta_text)

    else:
        fasta_text = st.sidebar.text_area(
            "Paste FASTA sequences:",
            height=200,
            placeholder=">sequence1\nATGGCTAGC...\n>sequence2\nATGTTAGCC..."
        )

        if fasta_text.strip():
            headers, sequences = parse_fasta_from_text(fasta_text)

    if sequences and len(sequences) > 0:
        st.sidebar.success(f"✅ {len(sequences)} sequences loaded")

        # Analysis options
        st.sidebar.header("🔬 Analysis Options")
        analysis_type = st.sidebar.selectbox(
            "Select analysis type:",
            [
                "Transcriptome-wide RSCU",
                "Per-gene RSCU",
                "Amino Acid Usage",
                "Codon Usage per 1000",
                "Relative Codon Frequencies"
            ]
        )

        if st.sidebar.button("🚀 Run Analysis"):
            with st.spinner("Analyzing sequences..."):

                if analysis_type == "Transcriptome-wide RSCU":
                    st.header("📊 Transcriptome-wide RSCU Analysis")

                    # Compute codon frequencies and RSCU
                    df_codcount, skipped = compute_codon_frequencies(headers, sequences)
                    df_rscu = compute_rscu_weights(df_codcount)

                    if skipped:
                        st.warning(f"⚠️ Skipped {len(skipped)} sequences (not multiple of 3): {', '.join(skipped[:5])}")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.subheader("📈 Codon Usage Frequencies")
                        fig_bar = create_codon_usage_bar_plot(df_codcount)
                        st.plotly_chart(fig_bar, use_container_width=True)

                    with col2:
                        st.subheader("📊 RSCU Distribution")
                        fig_hist = create_rscu_distribution(df_rscu)
                        st.plotly_chart(fig_hist, use_container_width=True)

                    st.subheader("🔥 RSCU Heatmap")
                    fig_heatmap = create_rscu_heatmap(df_rscu)
                    st.plotly_chart(fig_heatmap, use_container_width=True)

                    # Display data tables
                    st.subheader("📋 Detailed Results")

                    tab1, tab2 = st.tabs(["RSCU Values", "Codon Frequencies"])

                    with tab1:
                        st.dataframe(df_rscu, use_container_width=True)

                        # Download button
                        csv_rscu = df_rscu.to_csv(index=False)
                        st.download_button(
                            label="📥 Download RSCU Results (CSV)",
                            data=csv_rscu,
                            file_name="rscu_results.csv",
                            mime="text/csv"
                        )

                    with tab2:
                        st.dataframe(df_codcount, use_container_width=True)

                        csv_freq = df_codcount.to_csv(index=False)
                        st.download_button(
                            label="📥 Download Codon Frequencies (CSV)",
                            data=csv_freq,
                            file_name="codon_frequencies.csv",
                            mime="text/csv"
                        )

                elif analysis_type == "Per-gene RSCU":
                    st.header("🧬 Per-gene RSCU Analysis")

                    df_per_gene, skipped = compute_per_gene_rscu(headers, sequences)

                    if skipped:
                        st.warning(f"⚠️ Skipped {len(skipped)} sequences: {', '.join(skipped[:5])}")

                    if not df_per_gene.empty:
                        st.subheader("📊 RSCU Summary Statistics")

                        # Summary statistics
                        summary_stats = df_per_gene.groupby('Codon')['RSCU'].agg([
                            'count', 'mean', 'std', 'min', 'max'
                        ]).round(3)
                        st.dataframe(summary_stats, use_container_width=True)

                        st.subheader("📋 Per-gene RSCU Values")
                        st.dataframe(df_per_gene, use_container_width=True)

                        csv_per_gene = df_per_gene.to_csv(index=False)
                        st.download_button(
                            label="📥 Download Per-gene RSCU (CSV)",
                            data=csv_per_gene,
                            file_name="per_gene_rscu.csv",
                            mime="text/csv"
                        )

                elif analysis_type == "Amino Acid Usage":
                    st.header("🔬 Amino Acid Usage Analysis")

                    df_codcount, skipped = compute_codon_frequencies(headers, sequences)
                    df_rscu = compute_rscu_weights(df_codcount)
                    aa_df = compute_amino_acid_usage(df_rscu)

                    if skipped:
                        st.warning(f"⚠️ Skipped {len(skipped)} sequences: {', '.join(skipped[:5])}")

                    col1, col2 = st.columns(2)

                    with col1:
                        st.subheader("📊 Expected vs Observed Usage")
                        fig_aa = create_aa_usage_comparison(aa_df)
                        st.plotly_chart(fig_aa, use_container_width=True)

                    with col2:
                        st.subheader("📈 Usage Deviation")
                        aa_df['Deviation'] = aa_df['Obs_Freq(%)'] - aa_df['Expected_Freq(%)']
                        fig_dev = px.bar(
                            aa_df,
                            x='Amino_acid',
                            y='Deviation',
                            title="Deviation from Expected Usage",
                            color='Deviation',
                            color_continuous_scale='RdBu'
                        )
                        st.plotly_chart(fig_dev, use_container_width=True)

                    st.subheader("📋 Amino Acid Usage Results")
                    st.dataframe(aa_df, use_container_width=True)

                    csv_aa = aa_df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Amino Acid Usage (CSV)",
                        data=csv_aa,
                        file_name="amino_acid_usage.csv",
                        mime="text/csv"
                    )

                elif analysis_type == "Codon Usage per 1000":
                    st.header("📊 Codon Usage per 1000")

                    df_cu1000, skipped = compute_codon_usage_per_1000(headers, sequences)

                    if skipped:
                        st.warning(f"⚠️ Skipped {len(skipped)} sequences: {', '.join(skipped[:5])}")

                    fig_cu1000 = px.bar(
                        df_cu1000[df_cu1000['Amino_Acid'] != 'STOP'],
                        x='Codon',
                        y='Usage_per_1000',
                        color='Amino_Acid',
                        title="Codon Usage per 1000 Codons"
                    )
                    fig_cu1000.update_layout(height=500)
                    st.plotly_chart(fig_cu1000, use_container_width=True)

                    st.subheader("📋 Results")
                    st.dataframe(df_cu1000, use_container_width=True)

                    csv_cu1000 = df_cu1000.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Codon Usage per 1000 (CSV)",
                        data=csv_cu1000,
                        file_name="codon_usage_per_1000.csv",
                        mime="text/csv"
                    )

                elif analysis_type == "Relative Codon Frequencies":
                    st.header("📊 Relative Codon Frequencies")

                    df_rel_freq, skipped = compute_relative_codon_frequencies(headers, sequences)

                    if skipped:
                        st.warning(f"⚠️ Skipped {len(skipped)} sequences: {', '.join(skipped[:5])}")

                    if not df_rel_freq.empty:
                        st.subheader("📋 Relative Codon Frequencies per Gene")
                        st.dataframe(df_rel_freq, use_container_width=True)

                        csv_rel_freq = df_rel_freq.to_csv(index=False)
                        st.download_button(
                            label="📥 Download Relative Frequencies (CSV)",
                            data=csv_rel_freq,
                            file_name="relative_codon_frequencies.csv",
                            mime="text/csv"
                        )

    else:
        st.info("👆 Please upload a FASTA file or paste FASTA sequences in the sidebar to begin analysis.")

        # Show example
        st.markdown("### 📝 Example FASTA format:")
        st.code("""
>gene1
ATGGCTAAGTAG
>gene2
ATGTTTGCCTAG
>gene3
ATGGCCAAATAG
        """, language="text")

        st.markdown("""
        ### 🔬 Available Analysis Types:

        - **Transcriptome-wide RSCU**: Analyze codon usage across all sequences
        - **Per-gene RSCU**: Calculate RSCU for each individual gene
        - **Amino Acid Usage**: Compare expected vs observed amino acid frequencies
        - **Codon Usage per 1000**: Normalize codon usage per 1000 codons
        - **Relative Codon Frequencies**: Calculate relative frequencies for each gene
        """)


if __name__ == "__main__":
    main()