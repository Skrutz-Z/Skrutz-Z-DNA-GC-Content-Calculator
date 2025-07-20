# 🧬 DNA Sequence GC Content Calculator & Fragment Analyzer

A comprehensive web application for analyzing DNA sequences, calculating GC content percentages, and performing automatic sequence fragmentation with Z-DNA propensity analysis. Built with Streamlit and designed for bioinformatics research and education.

**🌐 Live Application:** https://skrutz-z-dna-gc-content-calculator.streamlit.app/

## 🔐 Access
**Password Protected Application** - Contact the administrator for access credentials.

## ✨ Features

### 🔬 Core Analysis
- **GC Content Calculation**: Accurate percentage calculation of Guanine-Cytosine content
- **Z-DNA Propensity Analysis**: Advanced algorithm for predicting Z-DNA formation likelihood
- **Nucleotide Composition**: Detailed breakdown of A, T, G, C percentages and counts
- **Sequence Validation**: Automatic sanitization and validation of DNA sequences
- **Batch Processing**: Handle up to 1000 sequences simultaneously

### 🧩 Automatic Fragmentation System
- **200-Base Fragmentation**: Automatically splits sequences into 200-nucleotide chunks
- **Position Tracking**: Precise start and end position tracking for each fragment
- **Fragment Analysis**: Complete GC content and Z-DNA analysis for each fragment
- **Smart Naming**: Descriptive fragment names with position ranges (e.g., Gene1_F1[1-200])
- **Batch Fragment Processing**: Process multiple sequences with automatic fragmentation

### 📊 Multiple Input Formats
- **CSV/Excel Files**: Upload spreadsheets with 'Gene Name' and 'Sequence' columns
- **FASTA Files**: Direct support for standard FASTA format
- **GenBank Files**: Parse GenBank format with annotations
- **Manual Entry**: Type sequences directly with custom gene names
- **Random Sequence Generator**: Generate test sequences with file upload capability

### 📈 Interactive Visualizations
- **GC Content Distribution**: Bar charts showing GC content across sequences
- **Fragment Statistics Dashboard**: Key metrics for fragment analysis
- **Nucleotide Composition**: Detailed breakdown for each sequence and fragment
- **Interactive Plots**: Plotly-based interactive visualizations
- **Heatmaps**: Visual representation of nucleotide patterns
- **Sequence Logos**: Multiple sequence alignment visualization

### 💾 Enhanced Export Options
- **Original Sequences**: Export in Excel, FASTA, GenBank formats
- **Fragment Analysis**: Separate exports for all 200-base fragments
- **Combined Exports**: Original sequences + fragments in single files
- **Z-DNA Data**: All exports include Z-DNA propensity calculations
- **Position Information**: Fragment start/end positions preserved
- **ZIP Archives**: Complete dataset exports with all formats

### 🎲 Additional Tools
- **Random Sequence Generator**: Generate test sequences with file upload support
- **Sequence Fragmenter Tool**: Custom fragment length selection (25, 50, 100, 200, 400 bases)
- **Bioinformatics Glossary**: Educational resource with key terms
- **Fun Facts**: DNA and genomics trivia
- **Tutorial Mode**: Step-by-step usage guide
- **Password Protection**: Secure access control

## 🚀 Quick Start

### Local Installation
1. **Clone the repository:**
   ```bash
   git clone https://github.com/yourusername/dna-gc-calculator.git
   cd dna-gc-calculator
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the application:**
   ```bash
   streamlit run app.py
   ```

4. **Open your browser** and navigate to `http://localhost:8501`

### Online Deployment
The application is deployed on Streamlit Cloud and can be accessed at: [Your Streamlit Cloud URL]

## 📋 Requirements
- Python 3.7+
- Streamlit
- Pandas
- Matplotlib
- Biopython
- Plotly
- Seaborn
- OpenPyXL
- NumPy

## 🧪 Sample Data
The repository includes sample files for testing:
- `random_sequences.fasta`: Sample FASTA sequences
- `Nucleotide content - spereadsheet example.xlsx`: Sample Excel file with sequences

## 🔧 Usage Instructions

### 1. Authentication
- Enter the required password to access the application
- Use the logout button in the sidebar to reset authentication

### 2. Choose Input Method
- **Upload CSV/Excel**: Upload files with 'Gene Name' and 'Sequence' columns
- **Upload FASTA**: Direct FASTA file support
- **Upload GenBank**: GenBank format parsing
- **Manual Entry**: Direct sequence input
- **Random Generator**: Generate or upload test sequences

### 3. Analyze Sequences
- Click **"Calculate GC Content & Fragment Analysis"** to process sequences
- The application automatically:
  - Analyzes original sequences
  - Splits into 200-base fragments
  - Calculates GC content and Z-DNA propensity for each fragment
  - Displays comprehensive results

### 4. Explore Results
- **Original Sequences Table**: Complete analysis of input sequences
- **Fragment Analysis Table**: Detailed analysis of all 200-base fragments
- **Fragment Statistics**: Key metrics and summary statistics
- **Interactive Visualizations**: Charts and plots for data exploration

### 5. Export Results
- **Original Sequences**: Download in Excel, FASTA, or GenBank format
- **Fragment Analysis**: Separate downloads for fragment data
- **Combined Export**: All data in single Excel file
- **ZIP Archive**: Complete dataset with all formats

## 🧬 Bioinformatics Features

### GC Content Analysis
- Calculates GC content percentage for each sequence and fragment
- Identifies high/low GC content regions
- Provides statistical analysis of nucleotide distribution
- Tracks GC content variation across fragments

### Z-DNA Propensity Analysis
- **Advanced Algorithm**: Analyzes alternating purine-pyrimidine patterns
- **6-Base Window Analysis**: Sliding window approach for pattern detection
- **GC-Rich Region Detection**: Higher weighting for GC-rich alternating regions
- **Percentage Scoring**: Normalized propensity scores for easy interpretation

### Sequence Processing
- Automatic sequence sanitization (removes invalid characters)
- Case-insensitive input handling
- Support for various sequence formats
- Fragment position tracking and naming

### Fragment Analysis
- **Automatic 200-Base Splitting**: Consistent fragment generation
- **Position Tracking**: Start and end positions for each fragment
- **Parent Gene Association**: Links fragments to original sequences
- **Comprehensive Metrics**: GC content, Z-DNA, nucleotide counts for each fragment

### Educational Resources
- Comprehensive glossary of bioinformatics terms
- Interactive tutorials and guides
- DNA and genomics fun facts
- Step-by-step usage instructions

## 🎯 Key Enhancements

### Automatic Fragmentation
- **Default 200-base chunks** for consistent analysis
- **Smart fragment naming** with position ranges
- **Complete analysis** for each fragment
- **Batch processing** for multiple sequences

### Z-DNA Analysis
- **Novel algorithm** for Z-DNA propensity prediction
- **Pattern recognition** in alternating sequences
- **GC-content weighting** for accuracy
- **Percentage output** for easy interpretation

### Enhanced Export System
- **Multiple format support** for all data types
- **Fragment-specific exports** with complete analysis
- **Combined data exports** for comprehensive analysis
- **Position information** preserved in all formats

### Improved User Interface
- **Overlay sidebar** for better content visibility
- **Professional styling** with dark blue theme
- **Fragment statistics dashboard** with key metrics
- **Organized export options** for different use cases

## 👨‍💻 Developer Information
- **Created by**: Shubh Rakesh Nahar
- **Institution**: Troy University
- **Purpose**: Educational and research tool for DNA sequence analysis
- **Version**: Enhanced with fragmentation and Z-DNA analysis

## 📄 License
This project is developed for educational and research purposes at Troy University.

## 🤝 Contributing
This is an academic project. For questions or suggestions, please contact the developer.

## 📞 Support
For technical support or access requests, please contact the developer.

---

**🔬 Advanced Features Summary:**
- ✅ Automatic 200-base sequence fragmentation
- ✅ Z-DNA propensity calculation algorithm
- ✅ Comprehensive fragment analysis with position tracking
- ✅ Enhanced export options for all data types
- ✅ Professional overlay sidebar interface
- ✅ Password protection and secure access
- ✅ Batch processing for multiple sequences
- ✅ Educational resources and tutorials 
