import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from Bio import SeqIO
import json
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import numpy as np
from matplotlib.patches import Rectangle
import math
import zipfile
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io
import openpyxl

# --- Instructions ---
def show_instructions():
    st.markdown("""
    ### 📋 Instructions
    
    1. **Choose Input Method**:
       - **Upload CSV/Excel**: Upload a CSV or Excel file with 'Gene Name' and 'Sequence' columns
       - **Upload FASTA**: Upload a FASTA file containing DNA sequences
       - **Upload GenBank**: Upload a GenBank file containing DNA sequences
       - **Manual Entry**: Enter sequences manually with gene names
    
    2. **Input Format**:
       - Only A, T, G, C nucleotides are accepted (case insensitive)
       - Invalid characters will be automatically removed
       - Each sequence must have a unique gene name
    
    3. **Analysis**:
       - Click "Calculate GC Content" to process sequences
       - View results in the interactive table
       - Explore visualizations of GC content and nucleotide composition
    
    4. **Export Options**:
       - Export results in multiple formats (Excel, CSV, JSON)
       - Use "Export All Results" for a complete dataset export
       - Customize output filename before downloading
    """)

# --- Sequence Sanitization and Validation ---
def sanitize_sequence(seq):
    # Convert to string if not already
    if not isinstance(seq, str):
        seq = str(seq)
    return ''.join(filter(lambda x: x.upper() in ['A', 'T', 'G', 'C'], seq.upper()))

def is_valid_sequence(seq):
    # Convert to string if not already
    if not isinstance(seq, str):
        seq = str(seq)
    return all(base in ['A', 'T', 'G', 'C'] for base in seq.upper())

# --- Nucleotide Analysis ---
def analyze_sequence(name, seq):
    # Convert sequence to string if it's not already
    if not isinstance(seq, str):
        seq = str(seq)
    
    seq = sanitize_sequence(seq)
    length = len(seq)
    if length == 0:
        return {"Gene Name": name, "Error": "Invalid sequence or empty after sanitization"}

    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')
    gc_count = g_count + c_count
    at_count = a_count + t_count

    return {
        "Gene Name": name,
        "Sequence": seq,
        "Length": length,
        "A Count": a_count,
        "T Count": t_count,
        "G Count": g_count,
        "C Count": c_count,
        "A %": round((a_count / length) * 100, 2),
        "T %": round((t_count / length) * 100, 2),
        "G %": round((g_count / length) * 100, 2),
        "C %": round((c_count / length) * 100, 2),
        "GC %": round((gc_count / length) * 100, 2),
        "AT %": round((at_count / length) * 100, 2),
    }

def process_fasta(file, max_sequences=1000):
    sequences = []
    try:
        # Ensure the file is read from the beginning
        file.seek(0)
        # Read binary and decode to text for FASTA
        text = file.read().decode("utf-8")
        text_io = io.StringIO(text)
        records = SeqIO.parse(text_io, "fasta")
        count = 0
        for record in records:
            if count >= max_sequences:
                break
            # Convert Seq object to string
            sequences.append((str(record.id), str(record.seq)))
            count += 1
        if count == 0:
            st.warning("No sequences found in the uploaded FASTA file after parsing.")
        elif count >= max_sequences:
             st.info(f"Processed the first {max_sequences} sequences from the FASTA file.")
    except Exception as e:
        st.error(f"Error parsing FASTA file: {e}")
        # Return empty list on error
        return []

    return sequences

def process_genbank(file, max_sequences=1000):
    sequences = []
    try:
        # Ensure the file is read from the beginning
        file.seek(0)
        # Read binary and decode to text for GenBank
        text = file.read().decode("utf-8")
        text_io = io.StringIO(text)
        records = SeqIO.parse(text_io, "genbank")
        count = 0
        for record in records:
             if count >= max_sequences:
                 break
             # Ensure molecule_type is set for GenBank export compatibility later
             # For uploaded GenBank, try to get molecule_type from annotations, default to DNA
             mol_type = record.annotations.get("molecule_type", "DNA")
             sequences.append((str(record.id), str(record.seq)))
             count += 1
        if count == 0:
             st.warning("No sequences found in the uploaded GenBank file after parsing.")
        elif count >= max_sequences:
             st.info(f"Processed the first {max_sequences} sequences from the GenBank file.")

    except Exception as e:
        st.error(f"Error parsing GenBank file: {e}")
        # Return empty list on error
        return []

    return sequences

def display_visuals(df):
    st.subheader("📊 GC Content Distribution")
    fig, ax = plt.subplots()
    df.plot.bar(x='Gene Name', y='GC %', ax=ax, color='purple', legend=False)
    plt.ylabel('GC %')
    plt.xticks(rotation=45, ha='right')
    st.pyplot(fig)

    st.subheader("🧬 Nucleotide Composition Per Sequence")
    for _, row in df.iterrows():
        fig, ax = plt.subplots()
        ax.pie(
            [row['A %'], row['T %'], row['G %'], row['C %']],
            labels=['A %', 'T %', 'G %', 'C %'],
            autopct='%1.1f%%',
            startangle=90
        )
        ax.set_title(f"{row['Gene Name']} - Base % Composition")
        st.pyplot(fig)

def export_data(df, format_type):
    signature = "Made by Shubh Rakesh Nahar / Troy University"
    if format_type == "Excel":
        towrite = BytesIO()
        with pd.ExcelWriter(towrite, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='GC Content Analysis')
            # Get the workbook and the worksheet
            workbook = writer.book
            worksheet = writer.sheets['GC Content Analysis']
            
            # Add signature as a comment in cell A1
            worksheet.cell(row=1, column=1).comment = openpyxl.comments.Comment(signature, "Author")
            
            # Auto-adjust columns' width
            for column in worksheet.columns:
                max_length = 0
                column = [cell for cell in column]
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = (max_length + 2)
                worksheet.column_dimensions[column[0].column_letter].width = adjusted_width
            
            # Add some basic formatting
            for row in worksheet.iter_rows(min_row=1, max_row=1):
                for cell in row:
                    cell.font = cell.font.copy(bold=True)
            
            # Add signature in the last row
            last_row = len(df) + 2
            worksheet.cell(row=last_row, column=1, value=signature)
        
        towrite.seek(0)
        return towrite, "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx"
    elif format_type == "CSV":
        towrite = BytesIO()
        # Add signature as a comment in the first line
        towrite.write(f"# {signature}\n".encode())
        df.to_csv(towrite, index=False)
        towrite.seek(0)
        return towrite, "text/csv", "csv"
    else:  # JSON
        towrite = BytesIO()
        # Add signature to the JSON metadata
        json_data = {
            "metadata": {
                "signature": signature,
                "timestamp": pd.Timestamp.now().isoformat()
            },
            "data": df.to_dict(orient='records')
        }
        json_str = json.dumps(json_data, indent=2)
        towrite.write(json_str.encode())
        towrite.seek(0)
        return towrite, "application/json", "json"

def create_gc_heatmap(df):
    st.subheader("🌡️ GC Content Heatmap")
    
    # Create a matrix of GC content for each position
    sequences = df['Sequence'].tolist()
    max_len = max(len(seq) for seq in sequences)
    
    # Initialize matrix
    gc_matrix = np.zeros((len(sequences), max_len))
    
    # Fill matrix with GC content for each position
    for i, seq in enumerate(sequences):
        for j in range(len(seq)):
            if j < len(seq):
                window = seq[max(0, j-10):min(len(seq), j+11)]
                gc_count = window.count('G') + window.count('C')
                gc_matrix[i, j] = (gc_count / len(window)) * 100
    
    # Create heatmap using plotly
    fig = go.Figure(data=go.Heatmap(
        z=gc_matrix,
        x=list(range(max_len)),
        y=df['Gene Name'].tolist(),
        colorscale='Viridis',
        colorbar=dict(title='GC %')
    ))
    
    fig.update_layout(
        title='GC Content Distribution Across Sequences',
        xaxis_title='Position',
        yaxis_title='Gene Name',
        height=400 + (len(sequences) * 20)  # Adjust height based on number of sequences
    )
    
    st.plotly_chart(fig, use_container_width=True)

def calculate_information_content(freq):
    """Calculate information content in bits."""
    if freq == 0:
        return 0
    return freq * math.log2(freq * 4)  # 4 for number of nucleotides

def create_sequence_logo(sequences):
    st.subheader("🎨 Sequence Logo")
    
    # Calculate position frequency matrix
    max_len = max(len(seq) for seq in sequences)
    pfm = np.zeros((4, max_len))  # 4 nucleotides
    
    for seq in sequences:
        for i, base in enumerate(seq):
            if i < max_len:
                if base == 'A':
                    pfm[0, i] += 1
                elif base == 'T':
                    pfm[1, i] += 1
                elif base == 'G':
                    pfm[2, i] += 1
                elif base == 'C':
                    pfm[3, i] += 1
    
    # Normalize
    pfm = pfm / len(sequences)
    
    # Calculate information content
    ic = np.zeros(max_len)
    for i in range(max_len):
        ic[i] = sum(calculate_information_content(freq) for freq in pfm[:, i])
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Colors for nucleotides
    colors = {'A': '#2ecc71', 'T': '#e74c3c', 'G': '#f1c40f', 'C': '#3498db'}
    bases = ['A', 'T', 'G', 'C']
    
    # Plot each position
    for i in range(max_len):
        # Sort frequencies for this position
        freqs = pfm[:, i]
        sorted_indices = np.argsort(freqs)
        
        # Plot each base
        y_bottom = 0
        for idx in sorted_indices:
            if freqs[idx] > 0:
                height = freqs[idx] * ic[i]
                rect = Rectangle((i, y_bottom), 1, height,
                               facecolor=colors[bases[idx]],
                               edgecolor='black',
                               linewidth=0.5)
                ax.add_patch(rect)
                y_bottom += height
    
    # Customize the plot
    ax.set_xlim(0, max_len)
    ax.set_ylim(0, max(ic) * 1.1)
    ax.set_xlabel('Position')
    ax.set_ylabel('Bits')
    ax.set_title('Sequence Logo')
    
    # Add legend
    legend_elements = [Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black')
                      for color in colors.values()]
    ax.legend(legend_elements, bases, loc='upper right')
    
    # Remove spines
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    # Show the plot
    st.pyplot(fig)
    plt.close(fig)

def create_interactive_plots(df):
    st.subheader("📊 Interactive Plots")
    
    # GC Content Distribution
    fig_gc = px.box(df, y='GC %', title='GC Content Distribution',
                    hover_data=['Gene Name', 'Length', 'GC %'])
    fig_gc.update_traces(marker_color='purple')
    st.plotly_chart(fig_gc, use_container_width=True)
    
    # Nucleotide Composition
    fig_comp = px.bar(df, 
                     x='Gene Name',
                     y=['A %', 'T %', 'G %', 'C %'],
                     title='Nucleotide Composition by Gene',
                     barmode='group',
                     hover_data=['Length', 'GC %'])
    st.plotly_chart(fig_comp, use_container_width=True)
    
    # Length vs GC Content
    fig_scatter = px.scatter(df,
                           x='Length',
                           y='GC %',
                           color='GC %',
                           hover_data=['Gene Name', 'Length', 'GC %'],
                           title='Sequence Length vs GC Content')
    st.plotly_chart(fig_scatter, use_container_width=True)

def export_all_to_zip(df, output_filename):
    signature = "Made by Shubh Rakesh Nahar / Troy University"
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zip_file:
        # Excel with formatting
        excel_buffer = BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
            df.to_excel(writer, index=False, sheet_name='GC Content Analysis')
            # Get the workbook and the worksheet
            workbook = writer.book
            worksheet = writer.sheets['GC Content Analysis']
            
            # Add signature as a comment in cell A1
            worksheet.cell(row=1, column=1).comment = openpyxl.comments.Comment(signature, "Author")
            
            # Auto-adjust columns' width
            for column in worksheet.columns:
                max_length = 0
                column = [cell for cell in column]
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = (max_length + 2)
                worksheet.column_dimensions[column[0].column_letter].width = adjusted_width
            
            # Add some basic formatting
            for row in worksheet.iter_rows(min_row=1, max_row=1):
                for cell in row:
                    cell.font = cell.font.copy(bold=True)
            
            # Add signature in the last row
            last_row = len(df) + 2
            worksheet.cell(row=last_row, column=1, value=signature)
        
        excel_buffer.seek(0)
        zip_file.writestr(f"{output_filename}.xlsx", excel_buffer.read())
        
        # CSV
        csv_buffer = BytesIO()
        # Add signature as a comment in the first line
        csv_buffer.write(f"# {signature}\n".encode())
        df.to_csv(csv_buffer, index=False)
        csv_buffer.seek(0)
        zip_file.writestr(f"{output_filename}.csv", csv_buffer.read())
        
        # JSON
        json_data = {
            "metadata": {
                "signature": signature,
                "timestamp": pd.Timestamp.now().isoformat()
            },
            "data": df.to_dict(orient='records')
        }
        json_str = json.dumps(json_data, indent=2)
        json_buffer = BytesIO()
        json_buffer.write(json_str.encode())
        json_buffer.seek(0)
        zip_file.writestr(f"{output_filename}.json", json_buffer.read())
    
    zip_buffer.seek(0)
    return zip_buffer

def export_to_fasta(sequences, result_df):
    """Export sequences to FASTA format with analysis results in description"""
    signature = "Made by Shubh Rakesh Nahar / Troy University"
    output = io.StringIO()
    # Add signature as a comment at the top
    output.write(f"# {signature}\n")
    for idx, row in result_df.iterrows():
        # Create a description with analysis results
        description = f"GC Content: {row['GC %']:.2f}%, Length: {row['Length']}, "
        description += f"A: {row['A Count']}, T: {row['T Count']}, G: {row['G Count']}, C: {row['C Count']}"
        
        # Write sequence with description
        output.write(f">{row['Gene Name']} {description}\n")
        output.write(f"{row['Sequence']}\n")
    
    # Add signature at the end
    output.write(f"\n# {signature}\n")
    return output.getvalue()

def export_to_genbank(sequences, result_df):
    """Export sequences to GenBank format with analysis results in features"""
    signature = "Made by Shubh Rakesh Nahar / Troy University"
    output = io.StringIO()
    for idx, row in result_df.iterrows():
        # Write GenBank header
        output.write(f"LOCUS       {row['Gene Name']}              {row['Length']} bp    DNA     linear\n")
        output.write(f"DEFINITION  {row['Gene Name']}\n")
        output.write(f"ACCESSION   {row['Gene Name']}\n")
        output.write(f"VERSION     {row['Gene Name']}\n")
        output.write(f"SOURCE      .\n")
        output.write(f"  ORGANISM  .\n")
        output.write(f"            .\n")
        
        # Add signature in the features section
        output.write("FEATURES             Location/Qualifiers\n")
        output.write("     source          1..{}\n".format(row['Length']))
        output.write('                     /organism="."\n')
        output.write('                     /mol_type="genomic DNA"\n')
        output.write('                     /note="{}"\n'.format(signature))
        output.write('     misc_feature    1..{}\n'.format(row['Length']))
        output.write('                     /note="GC Content: {:.2f}%"\n'.format(row['GC %']))
        output.write('                     /note="A Count: {}"\n'.format(row['A Count']))
        output.write('                     /note="T Count: {}"\n'.format(row['T Count']))
        output.write('                     /note="G Count: {}"\n'.format(row['G Count']))
        output.write('                     /note="C Count: {}"\n'.format(row['C Count']))
        
        # Write sequence
        output.write("ORIGIN\n")
        sequence = row['Sequence']
        for i in range(0, len(sequence), 60):
            chunk = sequence[i:i+60]
            output.write(f"     {i+1:9d} {chunk}\n")
        output.write("//\n")
    
    return output.getvalue()

def main():
    st.set_page_config(page_title="GC Content Calculator", page_icon=None, layout="wide")
    
    # Password protection
    if 'authenticated' not in st.session_state:
        st.session_state['authenticated'] = False
    
    if not st.session_state['authenticated']:
        st.markdown("""
            <h1 style='text-align: center; color: #6c63ff;'>🔐 DNA Sequence GC Content Calculator</h1>
            <h3 style='text-align: center; color: #ff6b6b;'>Authentication Required</h3>
        """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Create a centered login form
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            st.markdown("""
                <div style='text-align: center; padding: 20px; border: 2px solid #6c63ff; border-radius: 10px; background-color: #f8f9fa;'>
                    <h4>🔒 Please enter the password to access the application</h4>
                </div>
            """, unsafe_allow_html=True)
            
            password = st.text_input("Password", type="password", key="password_input")
            
            col_a, col_b, col_c = st.columns([1, 1, 1])
            with col_b:
                if st.button("🔓 Login", use_container_width=True):
                    if password == "TroyDNA2024":
                        st.session_state['authenticated'] = True
                        st.rerun()
                    else:
                        st.error("❌ Incorrect password. Please try again.")
            
            st.markdown("---")
            st.info("💡 **Hint**: The password is related to Troy University and DNA analysis.")
            
            # Add some DNA-themed styling
            st.markdown("""
                <div style='text-align: center; margin-top: 30px;'>
                    <p style='color: #6c63ff; font-style: italic;'>
                        🧬 Unlock the power of DNA sequence analysis 🧬
                    </p>
                </div>
            """, unsafe_allow_html=True)
        
        return
    
    # Main application content (only shown after authentication)
    st.markdown("""
        <h1 style='text-align: center; color: #6c63ff;'>🧬 GC Content Calculator</h1>
    """, unsafe_allow_html=True)
    
    # Add logout button in sidebar
    if st.sidebar.button("🚪 Logout"):
        st.session_state['authenticated'] = False
        st.rerun()

    sequences = []  # Ensure sequences is always defined

    # Sidebar with logo, navigation, and about info
    st.sidebar.image("https://cdn-icons-png.flaticon.com/512/616/616494.png", width=100)
    st.sidebar.title("Navigation")
    with st.sidebar.expander("About this app"):
        st.write("""
        This tool calculates GC content and provides interactive visualizations for DNA sequences. 
        Upload your data or enter sequences manually to explore nucleotide composition, GC content, and more. 
        
        **Made by Shubh Rakesh Nahar, Troy University.**
        """)
    # Fun facts about DNA/genes/sequences
    facts = [
        "The human genome contains about 3 billion base pairs.",
        "GC content can affect the stability of DNA.",
        "Some bacteria have extremely high or low GC content.",
        "DNA was first isolated by Friedrich Miescher in 1869.",
        "GC-rich regions are often found near gene promoters.",
        "Genes are segments of DNA that code for proteins.",
        "The longest human gene is over 2.4 million base pairs long!",
        "Mitochondrial DNA is inherited only from your mother.",
        "Some viruses use RNA instead of DNA as their genetic material.",
        "The fruit fly has about 15,000 genes, while humans have about 20,000-25,000.",
        "DNA stands for Deoxyribonucleic Acid.",
        "The double helix structure of DNA was discovered in 1953.",
        "Some plants have much more DNA than humans!",
        "The GC content of a genome can be used to identify species.",
        "DNA can be extracted from almost any living thing, even ancient fossils!"
    ]
    if 'fun_fact_idx' not in st.session_state:
        st.session_state['fun_fact_idx'] = random.randint(0, len(facts)-1)
    if st.sidebar.button("Show another fun fact"):
        prev_idx = st.session_state['fun_fact_idx']
        new_idx = prev_idx
        while new_idx == prev_idx:
            new_idx = random.randint(0, len(facts)-1)
        st.session_state['fun_fact_idx'] = new_idx
        st.rerun()
    st.sidebar.success(f"Fun Fact: {facts[st.session_state['fun_fact_idx']]}")

    # Glossary of Bioinformatics Terms
    glossary = {
        "GC Content": "The percentage of bases in a DNA or RNA molecule that are either guanine (G) or cytosine (C).",
        "FASTA": "A text-based format for representing nucleotide or peptide sequences.",
        "GenBank": "A rich file format for DNA sequences with annotations, used by NCBI.",
        "ORF": "Open Reading Frame, a sequence of DNA that could potentially encode a protein.",
        "Codon": "A sequence of three nucleotides that together form a unit of genetic code.",
        "Motif": "A short, recurring pattern in DNA that is presumed to have a biological function.",
        "CpG Island": "A region with a high frequency of CG dinucleotides, often found near gene promoters.",
        "SNP": "Single Nucleotide Polymorphism, a variation at a single position in a DNA sequence among individuals.",
        "BLAST": "A tool for comparing an input sequence against a database of sequences.",
        "Reverse Complement": "The sequence formed by reversing a DNA sequence and replacing each base with its complement.",
        "Translation": "The process of converting a nucleotide sequence into a protein sequence.",
        "Phylogenetic Tree": "A branching diagram showing evolutionary relationships among sequences.",
        "Alignment": "The arrangement of two or more sequences to identify regions of similarity."
    }
    with st.sidebar.expander("Glossary of Bioinformatics Terms"):
        search_term = st.text_input("Search glossary", key="glossary_search")
        for term, definition in glossary.items():
            if not search_term or search_term.lower() in term.lower():
                st.markdown(f"**{term}:** {definition}")

    # Tutorial Mode
    with st.sidebar.expander("Tutorial / Step-by-Step Guide"):
        st.markdown("""
        **How to use the GC Content Calculator:**
        1. **Choose Input Method:** Select how you want to input your sequences (CSV/Excel, FASTA, GenBank, or Manual Entry).
        2. **Upload or Enter Sequences:** Provide your DNA sequences using the chosen method.
        3. **Analyze:** Click 'Calculate GC Content' to process your sequences.
        4. **Explore Results:** View tables and interactive plots in the main area.
        5. **Export:** Download your results in Excel, CSV, JSON, FASTA, or GenBank format.
        6. **Fun Fact:** Click 'Show another fun fact' for a new DNA/genomics fact!
        """)

    # Random Sequence Generator
    with st.sidebar.expander("Random Sequence Generator"):
        st.markdown("Generate random DNA or protein sequences for practice or testing.")
        seq_type = st.selectbox("Sequence type", ["DNA", "Protein"], key="rand_seq_type")
        rand_num = st.number_input("Number of sequences", min_value=1, max_value=1000, value=1, key="rand_num")
        rand_len = st.number_input("Length of each sequence", min_value=5, max_value=40000, value=50, key="rand_len")
        if st.button("Generate Random Sequences"):
            import string
            import random as pyrandom
            if seq_type == "DNA":
                alphabet = "ATGC"
            else:
                alphabet = "ACDEFGHIKLMNPQRSTVWY"  # 20 amino acids
            
            # Generate sequences in chunks to manage memory
            chunk_size = 100  # Process 100 sequences at a time
            all_sequences = []
            
            for i in range(0, rand_num, chunk_size):
                current_chunk = min(chunk_size, rand_num - i)
                chunk_sequences = [(f"Random_{i+j+1}", ''.join(pyrandom.choices(alphabet, k=rand_len))) 
                                 for j in range(current_chunk)]
                all_sequences.extend(chunk_sequences)
                
                # Show progress
                progress = min(100, int((i + current_chunk) / rand_num * 100))
                st.progress(progress)
            
            if 'random_sequences' not in st.session_state:
                st.session_state['random_sequences'] = []
            st.session_state['random_sequences'] = all_sequences
            st.toast(f"{rand_num} random {seq_type} sequences of length {rand_len} generated!", icon=None)
        
        # Option to clear generated sequences
        if st.button("Clear Random Sequences"):
            st.session_state['random_sequences'] = []
            st.toast("Random sequences cleared.", icon=None)
        
        # Download random sequences if present
        if 'random_sequences' in st.session_state and st.session_state['random_sequences']:
            # First analyze the random sequences in chunks
            chunk_size = 100  # Process 100 sequences at a time
            all_results = []
            
            for i in range(0, len(st.session_state['random_sequences']), chunk_size):
                chunk = st.session_state['random_sequences'][i:i + chunk_size]
                chunk_results = [analyze_sequence(name, seq) for name, seq in chunk]
                all_results.extend(chunk_results)
                
                # Show progress
                progress = min(100, int((i + len(chunk)) / len(st.session_state['random_sequences']) * 100))
                st.progress(progress)
            
            rand_df = pd.DataFrame(all_results)
            
            # Excel with formatting
            excel_buffer = BytesIO()
            with pd.ExcelWriter(excel_buffer, engine='openpyxl') as writer:
                rand_df.to_excel(writer, index=False, sheet_name='GC Content Analysis')
                # Get the workbook and the worksheet
                workbook = writer.book
                worksheet = writer.sheets['GC Content Analysis']
                
                # Auto-adjust columns' width
                for column in worksheet.columns:
                    max_length = 0
                    column = [cell for cell in column]
                    for cell in column:
                        try:
                            if len(str(cell.value)) > max_length:
                                max_length = len(str(cell.value))
                        except:
                            pass
                    adjusted_width = (max_length + 2)
                    worksheet.column_dimensions[column[0].column_letter].width = adjusted_width
                
                # Add some basic formatting
                for row in worksheet.iter_rows(min_row=1, max_row=1):
                    for cell in row:
                        cell.font = cell.font.copy(bold=True)
            
            excel_buffer.seek(0)
            st.download_button(
                label="Download Random Sequences (Excel)",
                data=excel_buffer,
                file_name="random_sequences.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
            
            # FASTA
            fasta_buffer = export_to_fasta(st.session_state['random_sequences'], rand_df)
            st.download_button(
                label="Download Random Sequences (FASTA)",
                data=fasta_buffer,
                file_name="random_sequences.fasta",
                mime="text/plain"
            )
            
            # GenBank
            gb_buffer = export_to_genbank(st.session_state['random_sequences'], rand_df)
            st.download_button(
                label="Download Random Sequences (GenBank)",
                data=gb_buffer,
                file_name="random_sequences.gb",
                mime="text/plain"
            )

    # Add random sequences to main input if present
    if 'random_sequences' in st.session_state and st.session_state['random_sequences']:
        sequences.extend(st.session_state['random_sequences'])
        st.session_state['random_sequences'] = [] # Clear after adding

    st.markdown("---")
    st.markdown("#### Upload your data or enter sequences manually below.")
    st.info("You can input up to 1000 sequences at a time (CSV, Excel, FASTA, GenBank, or manual entry). If you upload more, only the first 1000 will be processed.")

    input_method = st.radio("Choose input method", ["Upload CSV/Excel", "Upload FASTA", "Upload GenBank", "Manual Entry"])

    if input_method == "Upload CSV/Excel":
        uploaded_file = st.file_uploader("Upload a CSV or Excel file with 'Gene Name' and 'Sequence' columns", type=["csv", "xlsx", "xls"])
        if uploaded_file:
            try:
                uploaded_file.seek(0) # Ensure pointer is at the start
                if uploaded_file.name.endswith((".xlsx", ".xls")):
                    df = pd.read_excel(uploaded_file)
                else:
                    df = pd.read_csv(uploaded_file)
                
                if "Gene Name" in df.columns and "Sequence" in df.columns:
                    # Convert sequences to strings and handle any NaN values
                    df['Sequence'] = df['Sequence'].fillna('').astype(str)
                    df['Gene Name'] = df['Gene Name'].fillna('').astype(str)
                    if len(df) > 1000:
                        st.toast("More than 1000 sequences detected. Only the first 1000 will be processed.", icon=None)
                    sequences = list(zip(df["Gene Name"], df["Sequence"]))[:1000]
                else:
                    st.toast("File must contain 'Gene Name' and 'Sequence' columns.", icon=None)
            except Exception as e:
                st.toast(f"Error reading file: {str(e)}", icon=None)

    elif input_method == "Upload FASTA":
        fasta_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])
        if fasta_file:
            fasta_file.seek(0)  # Reset pointer
            try:
                sequences = process_fasta(fasta_file, max_sequences=1000)
                if len(sequences) == 0:
                    st.warning("No sequences found in the uploaded FASTA file.")
                elif len(sequences) == 1000:
                    st.toast("More than 1000 sequences detected. Only the first 1000 will be processed.", icon=None)
            except Exception as e:
                st.toast(f"Error processing FASTA file: {str(e)}", icon=None)

    elif input_method == "Upload GenBank":
        gb_file = st.file_uploader("Upload a GenBank file", type=["gb", "gbk"])
        if gb_file:
            gb_file.seek(0)  # Reset pointer
            try:
                sequences = process_genbank(gb_file, max_sequences=1000)
                if len(sequences) == 0:
                    st.warning("No sequences found in the uploaded GenBank file.")
                elif len(sequences) == 1000:
                    st.toast("More than 1000 sequences detected. Only the first 1000 will be processed.", icon=None)
            except Exception as e:
                st.toast(f"Error processing GenBank file: {str(e)}", icon=None)

    elif input_method == "Manual Entry":
        num = st.number_input("How many sequences would you like to enter?", min_value=1, max_value=1000, value=1)
        manual_entries = []
        for i in range(num):
            st.markdown(f"**Sequence {i+1}**")
            name = st.text_input(f"Gene Name {i+1}", key=f"name_{i}")
            seq = st.text_area(f"Sequence {i+1}", key=f"seq_{i}")
            if name and seq:
                manual_entries.append((name, seq))
        sequences.extend(manual_entries)

    # Add random sequences to main input if present
    if 'random_sequences' in st.session_state and st.session_state['random_sequences']:
        sequences.extend(st.session_state['random_sequences'])
        st.session_state['random_sequences'] = [] # Clear after adding

    if sequences:
        if st.button("Calculate GC Content"):
            results = [analyze_sequence(name, seq) for name, seq in sequences]
            result_df = pd.DataFrame(results)

            if "Error" in result_df.columns:
                st.toast("Some sequences were invalid and skipped.", icon=None)
                result_df = result_df.dropna(subset=["Length"])

            # Check if there are valid sequences for output
            if len(result_df) == 0:
                 st.warning("No valid sequences were processed.")
            else:

                # Conditional Visualization Display
                if len(result_df) > 50:
                    st.info("Analysis complete. Visualizations and on-page data table are disabled for more than 50 sequences to ensure performance. Please download your full results below.")
                else:
                    st.toast("Analysis complete!", icon=None)
                    st.markdown("---")
                    st.markdown("#### Results & Visualizations")
                    
                    with st.expander("Show Data Table", expanded=True):
                        st.dataframe(result_df, use_container_width=True)
                    with st.expander("Show GC Content Heatmap"):
                        if len(result_df) > 0:
                            create_gc_heatmap(result_df)
                        else:
                            st.info("No data to plot.")
                    with st.expander("Show Sequence Logo"):
                        # Use the original 'sequences' list for the logo, not result_df
                        if len(sequences) > 1 and all(len(seq) > 0 for _, seq in sequences):
                            # Filter original sequences based on names in result_df to match analysis results
                            valid_gene_names = result_df['Gene Name'].tolist()
                            sequences_for_logo = [(name, seq) for name, seq in sequences if name in valid_gene_names]
                            if len(sequences_for_logo) > 1:
                                create_sequence_logo([seq for _, seq in sequences_for_logo])
                            else:
                                st.info("Not enough valid sequences to plot a sequence logo.")
                        else:
                            st.info("Not enough data to plot a sequence logo.")

                    with st.expander("Show Interactive Plots"):
                        if len(result_df) > 1:
                            create_interactive_plots(result_df)
                        else:
                            st.info("Not enough data to plot interactive plots.")

                # Export Results section - always shown if there are results
                st.markdown("This section should always appear if there are results.") # Added for debugging
                st.markdown("---") # Add a separator before export options
                st.subheader("Export Results")
                output_filename = st.text_input("Enter output file name (without extension):", "gc_output_v4", key="export_filename")
                
                # Export buttons
                if st.button("Export All Results"):
                    zip_buffer = export_all_to_zip(result_df, output_filename)
                    st.toast("All results exported as ZIP!", icon=None)
                    st.download_button(
                        label="Download All Results (ZIP)",
                        data=zip_buffer,
                        file_name=f"{output_filename}_all_results.zip",
                        mime="application/zip"
                    )
                excel_buffer = BytesIO()
                result_df.to_excel(excel_buffer, index=False, engine='openpyxl')
                excel_buffer.seek(0)
                st.download_button(
                    label="Download Excel Only",
                    data=excel_buffer,
                    file_name=f"{output_filename}.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )
                # FASTA export button
                # Use result_df for FASTA/GenBank exports to include only analyzed sequences
                fasta_sequences_to_export = [(row["Gene Name"], row["Sequence"]) for _, row in result_df.iterrows()]
                fasta_buffer = export_to_fasta(fasta_sequences_to_export, result_df)
                st.download_button(
                    label="Download as FASTA",
                    data=fasta_buffer,
                    file_name=f"{output_filename}.fasta",
                    mime="text/plain"
                )
                # GenBank export button
                gb_sequences_to_export = [(row["Gene Name"], row["Sequence"]) for _, row in result_df.iterrows()]
                gb_buffer = export_to_genbank(gb_sequences_to_export, result_df)
                st.download_button(
                    label="Download as GenBank",
                    data=gb_buffer,
                    file_name=f"{output_filename}.gb",
                    mime="text/plain"
                )

    st.markdown("---")
    st.markdown("<div style='text-align: center; color: #888;'>Made by Shubh Rakesh Nahar, Troy University </div>", unsafe_allow_html=True)

if __name__ == "__main__":
    main()
