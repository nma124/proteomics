import streamlit as st
import pandas as pd
import io

st.set_page_config(page_title="Peptide Data Analyzer", layout="wide")

st.title("Peptide Data Analyzer")
st.markdown("Upload your peptide data file (CSV or Excel) to view and download processed results")

# File upload
uploaded_file = st.file_uploader(
    "Choose a file", 
    type=['csv', 'xlsx', 'xls'],
    help="Upload CSV or Excel file containing peptide data"
)

if uploaded_file is not None:
    try:
        # Read file based on type
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        else:
            df = pd.read_excel(uploaded_file)
        
        # Display basic info
        st.success(f"✅ File loaded successfully: {uploaded_file.name}")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Rows", len(df))
        with col2:
            st.metric("Total Columns", len(df.columns))
        with col3:
            # Try to count unique peptides if column exists
            peptide_cols = [col for col in df.columns if 'peptide' in col.lower()]
            if peptide_cols:
                unique_peptides = df[peptide_cols[0]].nunique()
                st.metric("Unique Peptides", unique_peptides)
            else:
                st.metric("Unique Peptides", "N/A")
        
        # Display the data
        st.subheader("Data Preview")
        st.dataframe(df, use_container_width=True, height=400)
        
        # Basic statistics
        with st.expander("📊 Show Data Statistics"):
            st.write(df.describe())
        
        # Column info
        with st.expander("📋 Show Column Information"):
            col_info = pd.DataFrame({
                'Column': df.columns,
                'Type': df.dtypes,
                'Non-Null Count': df.count(),
                'Null Count': df.isnull().sum()
            })
            st.dataframe(col_info, use_container_width=True)
        
        # Download options
        st.subheader("Download Options")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Download as CSV
            csv = df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="📥 Download as CSV",
                data=csv,
                file_name=f"processed_{uploaded_file.name.rsplit('.', 1)[0]}.csv",
                mime="text/csv"
            )
        
        with col2:
            # Download as Excel
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine='xlsxwriter') as writer:
                df.to_excel(writer, index=False, sheet_name='Data')
            buffer.seek(0)
            
            st.download_button(
                label="📥 Download as Excel",
                data=buffer,
                file_name=f"processed_{uploaded_file.name.rsplit('.', 1)[0]}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        # Optional: Filter by peptide if peptide column exists
        if peptide_cols:
            st.subheader("Filter by Peptide")
            selected_peptides = st.multiselect(
                "Select peptides to display:",
                options=sorted(df[peptide_cols[0]].unique()),
                default=None
            )
            
            if selected_peptides:
                filtered_df = df[df[peptide_cols[0]].isin(selected_peptides)]
                st.write(f"Showing {len(filtered_df)} rows for {len(selected_peptides)} peptide(s)")
                st.dataframe(filtered_df, use_container_width=True)
                
                # Download filtered data
                csv_filtered = filtered_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📥 Download Filtered Data (CSV)",
                    data=csv_filtered,
                    file_name=f"filtered_{uploaded_file.name.rsplit('.', 1)[0]}.csv",
                    mime="text/csv",
                    key="filtered_download"
                )
    
    except Exception as e:
        st.error(f"❌ Error reading file: {str(e)}")
        st.info("Please ensure your file is properly formatted CSV or Excel file")

else:
    st.info("👆 Upload a file to get started")
    st.markdown("""
    ### Supported formats:
    - CSV (`.csv`)
    - Excel (`.xlsx`, `.xls`)
    
    ### Features:
    - View data in interactive table
    - Automatic peptide counting
    - Download as CSV or Excel
    - Filter by specific peptides
    - View statistics and column information
    """)
