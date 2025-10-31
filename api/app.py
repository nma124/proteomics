#!/usr/bin/env python3
"""
FastAPI Backend for Proteomics PRM Data Processing

Provides REST API endpoints for:
- File upload (MS data and concentration files)
- Data processing
- Results retrieval
- Job status tracking
"""

from fastapi import FastAPI, UploadFile, File, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse
from pydantic import BaseModel
from typing import Optional, Dict
import sys
from pathlib import Path
import uuid
import shutil
import json
from datetime import datetime
import pandas as pd
import numpy as np
import io
import base64
import os

# Matplotlib setup for headless rendering
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import linear_model

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))
from scripts.unified_processor import process_prm_unified

app = FastAPI(
    title="Proteomics PRM Processing API",
    description="Unified API for processing Skyline PRM data with automatic format detection",
    version="1.0.0"
)

# CORS middleware for React frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:5173"],  # React dev servers
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Storage for uploaded files and job status
UPLOAD_DIR = Path("api/uploads")
RESULTS_DIR = Path("api/results")
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# In-memory job tracking (use Redis/DB for production)
jobs: Dict[str, dict] = {}


class JobStatus(BaseModel):
    """Job status response model"""
    job_id: str
    status: str  # "pending", "processing", "completed", "failed"
    message: Optional[str] = None
    progress: Optional[int] = None
    result_file: Optional[str] = None
    created_at: str
    completed_at: Optional[str] = None


class ProcessRequest(BaseModel):
    """Request model for processing endpoint"""
    ms_file_id: str
    concentration_file_id: str
    output_filename: Optional[str] = "results.csv"


@app.get("/")
def root():
    """API health check"""
    return {
        "status": "online",
        "service": "Proteomics PRM Processing API",
        "version": "1.0.0",
        "endpoints": {
            "upload_ms": "/upload/ms",
            "upload_concentration": "/upload/concentration",
            "process": "/process",
            "job_status": "/jobs/{job_id}",
            "download": "/download/{job_id}"
        }
    }


@app.post("/upload/ms")
async def upload_ms_file(file: UploadFile = File(...)):
    """
    Upload mass spectrometry data file (Skyline export CSV)
    
    Returns:
        file_id: Unique identifier for the uploaded file
    """
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")
    
    # Generate unique file ID
    file_id = str(uuid.uuid4())
    file_path = UPLOAD_DIR / f"{file_id}_ms.csv"
    
    # Save uploaded file
    with file_path.open("wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    return {
        "file_id": file_id,
        "filename": file.filename,
        "size": file_path.stat().st_size,
        "type": "mass_spectrometry"
    }


@app.post("/upload/concentration")
async def upload_concentration_file(file: UploadFile = File(...)):
    """
    Upload peptide concentration/dilution CSV file
    
    Returns:
        file_id: Unique identifier for the uploaded file
    """
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")
    
    # Generate unique file ID
    file_id = str(uuid.uuid4())
    file_path = UPLOAD_DIR / f"{file_id}_conc.csv"
    
    # Save uploaded file
    with file_path.open("wb") as buffer:
        shutil.copyfileobj(file.file, buffer)
    
    return {
        "file_id": file_id,
        "filename": file.filename,
        "size": file_path.stat().st_size,
        "type": "concentration"
    }


def process_data_background(job_id: str, ms_file: str, conc_file: str, output_file: str):
    """Background task for data processing"""
    try:
        jobs[job_id]["status"] = "processing"
        jobs[job_id]["message"] = "Processing data..."
        jobs[job_id]["progress"] = 50
        
        # Run the unified processor
        result_df = process_prm_unified(ms_file, conc_file, output_file)
        
        jobs[job_id]["status"] = "completed"
        jobs[job_id]["message"] = f"Processing complete! Generated {result_df.shape[0]} rows × {result_df.shape[1]} columns"
        jobs[job_id]["progress"] = 100
        jobs[job_id]["result_file"] = output_file
        jobs[job_id]["completed_at"] = datetime.now().isoformat()
        jobs[job_id]["result_stats"] = {
            "rows": int(result_df.shape[0]),
            "columns": int(result_df.shape[1]),
            "peptides": int(result_df['Peptide'].nunique()) if 'Peptide' in result_df.columns else None
        }
        
    except Exception as e:
        jobs[job_id]["status"] = "failed"
        jobs[job_id]["message"] = f"Processing failed: {str(e)}"
        jobs[job_id]["progress"] = 0
        jobs[job_id]["completed_at"] = datetime.now().isoformat()


@app.post("/process", response_model=JobStatus)
async def process_data(request: ProcessRequest, background_tasks: BackgroundTasks):
    """
    Process uploaded data files
    
    Args:
        ms_file_id: File ID from MS file upload
        concentration_file_id: File ID from concentration file upload
        output_filename: Optional custom output filename
        
    Returns:
        job_id: Unique identifier for tracking job status
    """
    # Validate uploaded files exist
    ms_file = UPLOAD_DIR / f"{request.ms_file_id}_ms.csv"
    conc_file = UPLOAD_DIR / f"{request.concentration_file_id}_conc.csv"
    
    if not ms_file.exists():
        raise HTTPException(status_code=404, detail="MS file not found")
    if not conc_file.exists():
        raise HTTPException(status_code=404, detail="Concentration file not found")
    
    # Create job
    job_id = str(uuid.uuid4())
    output_file = RESULTS_DIR / f"{job_id}_{request.output_filename}"
    
    jobs[job_id] = {
        "job_id": job_id,
        "status": "pending",
        "message": "Job created, starting processing...",
        "progress": 0,
        "created_at": datetime.now().isoformat(),
        "ms_file": str(ms_file),
        "conc_file": str(conc_file),
        "output_file": str(output_file)
    }
    
    # Start background processing
    background_tasks.add_task(
        process_data_background,
        job_id,
        str(ms_file),
        str(conc_file),
        str(output_file)
    )
    
    return JobStatus(**jobs[job_id])


@app.get("/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(job_id: str):
    """
    Get job status and results
    
    Returns:
        Job status, progress, and result file path if completed
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    return JobStatus(**jobs[job_id])


@app.get("/download/{job_id}")
async def download_results(job_id: str):
    """
    Download processed results CSV
    
    Returns:
        CSV file with processed data
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    result_file = Path(job["result_file"])
    if not result_file.exists():
        raise HTTPException(status_code=404, detail="Result file not found")
    
    return FileResponse(
        path=result_file,
        filename=result_file.name,
        media_type="text/csv"
    )


@app.get("/jobs")
async def list_jobs():
    """List all jobs"""
    return {"jobs": list(jobs.values())}


@app.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    """Delete a job and its associated files"""
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    # Delete result file if exists
    if "result_file" in job and job["result_file"]:
        result_file = Path(job["result_file"])
        if result_file.exists():
            result_file.unlink()
    
    # Remove job from tracking
    del jobs[job_id]
    
    return {"message": "Job deleted successfully"}


@app.get("/analytics/{job_id}")
async def get_analytics(job_id: str):
    """
    Get analytics and statistics for processed data
    
    Returns regression statistics, quality metrics, and summary data
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    result_file = Path(job["result_file"])
    if not result_file.exists():
        raise HTTPException(status_code=404, detail="Result file not found")
    
    try:
        # Read the processed data
        df = pd.read_csv(result_file)
        df.columns = df.columns.str.strip().str.lower()
        
        # Extract regression statistics
        regression_cols = ['r2', 'intercept', 'gradient']
        agg_cols = ['mean_grad', 'stdv_grad', 'cov_grad', 'mean_intercept', 
                   'stdv_intercept', 'cov_intercept', 'mean_r2', 'stdv_r2', 'cov_r2']
        qtest_cols = ['qtest_grad', 'qtest_intercept', 'qtest_r2']
        
        # Check which columns are present
        present_regression = [col for col in regression_cols if col in df.columns]
        present_agg = [col for col in agg_cols if col in df.columns]
        present_qtest = [col for col in qtest_cols if col in df.columns]
        
        # Summary statistics
        summary = {
            "total_rows": int(df.shape[0]),
            "total_columns": int(df.shape[1]),
            "unique_peptides": int(df['peptide'].nunique()) if 'peptide' in df.columns else 0,
            "unique_fragments": int(df['fragment ion'].nunique()) if 'fragment ion' in df.columns else 0,
            "unique_plot_categories": int(df['plot_cat'].nunique()) if 'plot_cat' in df.columns else 0,
        }
        
        # Regression quality metrics (aggregate across all data)
        quality_metrics = {}
        if 'r2' in df.columns:
            r2_values = pd.to_numeric(df['r2'], errors='coerce').dropna()
            quality_metrics['r2_mean'] = float(r2_values.mean()) if len(r2_values) > 0 else None
            quality_metrics['r2_median'] = float(r2_values.median()) if len(r2_values) > 0 else None
            quality_metrics['r2_min'] = float(r2_values.min()) if len(r2_values) > 0 else None
            quality_metrics['r2_max'] = float(r2_values.max()) if len(r2_values) > 0 else None
            quality_metrics['r2_good_fits'] = int((r2_values > 0.95).sum()) if len(r2_values) > 0 else 0
            quality_metrics['r2_poor_fits'] = int((r2_values < 0.8).sum()) if len(r2_values) > 0 else 0
        
        # Per-peptide statistics
        peptide_stats = []
        if 'peptide' in df.columns and 'r2' in df.columns:
            for peptide, group in df.groupby('peptide'):
                r2_vals = pd.to_numeric(group['r2'], errors='coerce').dropna()
                peptide_stats.append({
                    "peptide": peptide,
                    "n_measurements": int(len(group)),
                    "n_fragments": int(group['fragment ion'].nunique()) if 'fragment ion' in group.columns else 0,
                    "r2_mean": float(r2_vals.mean()) if len(r2_vals) > 0 else None,
                    "r2_std": float(r2_vals.std()) if len(r2_vals) > 0 else None,
                    "gradient_mean": float(pd.to_numeric(group['gradient'], errors='coerce').mean()) if 'gradient' in group.columns else None,
                })
        
        # Get unique plot categories for plotting
        plot_categories = []
        if 'plot_cat' in df.columns:
            plot_categories = sorted(df['plot_cat'].dropna().unique().tolist())[:100]  # Limit to 100
        
        return {
            "job_id": job_id,
            "summary": summary,
            "quality_metrics": quality_metrics,
            "peptide_stats": peptide_stats,
            "plot_categories": plot_categories,
            "available_columns": {
                "regression": present_regression,
                "aggregates": present_agg,
                "qtest": present_qtest,
            }
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analytics computation failed: {str(e)}")


@app.get("/plots/{job_id}/{plot_category}")
async def generate_plot(job_id: str, plot_category: str):
    """
    Generate calibration plot for a specific plot category
    
    Returns base64-encoded PNG image
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    result_file = Path(job["result_file"])
    if not result_file.exists():
        raise HTTPException(status_code=404, detail="Result file not found")
    
    try:
        # Read data and filter for plot category
        df = pd.read_csv(result_file)
        df.columns = df.columns.str.strip().str.lower()
        
        if 'plot_cat' not in df.columns:
            raise HTTPException(status_code=400, detail="No plot categories in data")
        
        cat_df = df[df['plot_cat'] == plot_category]
        
        if cat_df.empty:
            raise HTTPException(status_code=404, detail=f"Plot category '{plot_category}' not found")
        
        # Clean data
        required_cols = ['heavy_conc', 'area_ratio']
        if not all(col in cat_df.columns for col in required_cols):
            raise HTTPException(status_code=400, detail="Missing required columns for plotting")
        
        cat_df = cat_df.replace([np.inf, -np.inf], np.nan).dropna(subset=required_cols)
        
        if cat_df.empty:
            raise HTTPException(status_code=400, detail="No valid data points for plotting")
        
        # Extract data
        x = cat_df['heavy_conc'].to_numpy(dtype=float)
        y = cat_df['area_ratio'].to_numpy(dtype=float)
        order = np.argsort(x)
        x_sorted, y_sorted = x[order], y[order]
        
        # Fit linear model
        TOL = 1e-6
        weights = np.where(np.abs(x_sorted) > TOL, 1.0 / x_sorted, 1.0)
        model = linear_model.LinearRegression()
        model.fit(x_sorted[:, np.newaxis], y_sorted[:, np.newaxis], sample_weight=weights)
        y_fit = model.predict(x_sorted[:, np.newaxis]).squeeze()
        intercept = float(model.intercept_.squeeze())
        slope = float(model.coef_.squeeze())
        
        # Calculate R²
        from sklearn.metrics import r2_score
        r2 = r2_score(y_sorted, y_fit)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(x_sorted, y_sorted, color='#2563eb', s=50, alpha=0.7, label='Observed')
        ax.plot(x_sorted, y_fit, color='#dc2626', linewidth=2, label='Linear Fit')
        ax.set_title(plot_category, fontsize=14, fontweight='bold')
        ax.set_xlabel('Heavy Concentration (ng/mL)', fontsize=11)
        ax.set_ylabel('Area Ratio (Light / Heavy)', fontsize=11)
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3, linestyle='--')
        
        # Add stats annotation
        stats_text = f"R² = {r2:.4f}\nSlope = {slope:.3g}\nIntercept = {intercept:.3g}"
        ax.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', 
                   va='top', ha='left', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        fig.tight_layout()
        
        # Convert to base64
        buffer = io.BytesIO()
        fig.savefig(buffer, format='png', dpi=100, bbox_inches='tight')
        plt.close(fig)
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.read()).decode('utf-8')
        
        return {
            "plot_category": plot_category,
            "image": f"data:image/png;base64,{img_base64}",
            "stats": {
                "r2": r2,
                "slope": slope,
                "intercept": intercept,
                "n_points": len(x_sorted)
            }
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Plot generation failed: {str(e)}")


@app.post("/plots/{job_id}/generate-all")
async def generate_all_plots(job_id: str, background_tasks: BackgroundTasks):
    """
    Generate all calibration plots for a job in the background
    
    Returns job status immediately and processes plots in background
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    result_file = Path(job["result_file"])
    if not result_file.exists():
        raise HTTPException(status_code=404, detail="Result file not found")
    
    # Create plots directory
    plots_dir = RESULTS_DIR / f"{job_id}_plots"
    plots_dir.mkdir(exist_ok=True)
    
    def generate_plots_background():
        try:
            df = pd.read_csv(result_file)
            df.columns = df.columns.str.strip().str.lower()
            
            if 'plot_cat' not in df.columns:
                return
            
            plot_categories = df['plot_cat'].dropna().unique()
            
            for idx, plot_cat in enumerate(plot_categories):
                try:
                    cat_df = df[df['plot_cat'] == plot_cat]
                    cat_df = cat_df.replace([np.inf, -np.inf], np.nan).dropna(subset=['heavy_conc', 'area_ratio'])
                    
                    if cat_df.empty:
                        continue
                    
                    x = cat_df['heavy_conc'].to_numpy(dtype=float)
                    y = cat_df['area_ratio'].to_numpy(dtype=float)
                    order = np.argsort(x)
                    x_sorted, y_sorted = x[order], y[order]
                    
                    # Fit model
                    TOL = 1e-6
                    weights = np.where(np.abs(x_sorted) > TOL, 1.0 / x_sorted, 1.0)
                    model = linear_model.LinearRegression()
                    model.fit(x_sorted[:, np.newaxis], y_sorted[:, np.newaxis], sample_weight=weights)
                    y_fit = model.predict(x_sorted[:, np.newaxis]).squeeze()
                    intercept = float(model.intercept_.squeeze())
                    slope = float(model.coef_.squeeze())
                    
                    from sklearn.metrics import r2_score
                    r2 = r2_score(y_sorted, y_fit)
                    
                    # Create plot
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.scatter(x_sorted, y_sorted, color='#2563eb', s=50, alpha=0.7, label='Observed')
                    ax.plot(x_sorted, y_fit, color='#dc2626', linewidth=2, label='Linear Fit')
                    ax.set_title(plot_cat, fontsize=14, fontweight='bold')
                    ax.set_xlabel('Heavy Concentration (ng/mL)', fontsize=11)
                    ax.set_ylabel('Area Ratio (Light / Heavy)', fontsize=11)
                    ax.legend(loc='best')
                    ax.grid(True, alpha=0.3, linestyle='--')
                    
                    stats_text = f"R² = {r2:.4f}\nSlope = {slope:.3g}\nIntercept = {intercept:.3g}"
                    ax.annotate(stats_text, xy=(0.05, 0.95), xycoords='axes fraction', 
                               va='top', ha='left', fontsize=10,
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
                    
                    fig.tight_layout()
                    
                    # Save plot
                    safe_name = plot_cat.replace('/', '-').replace(' ', '_')
                    plot_path = plots_dir / f"{safe_name}.png"
                    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
                    plt.close(fig)
                    
                except Exception as e:
                    print(f"Error generating plot for {plot_cat}: {e}")
                    continue
            
            # Mark plots as generated
            jobs[job_id]["plots_generated"] = True
            jobs[job_id]["plots_dir"] = str(plots_dir)
            jobs[job_id]["plots_count"] = len(list(plots_dir.glob("*.png")))
            
        except Exception as e:
            print(f"Error in bulk plot generation: {e}")
            jobs[job_id]["plots_generated"] = False
    
    # Start background task
    background_tasks.add_task(generate_plots_background)
    
    return {
        "job_id": job_id,
        "message": "Plot generation started in background",
        "status": "generating"
    }


@app.get("/plots/{job_id}/status")
async def get_plots_status(job_id: str):
    """
    Check if plots have been generated for a job
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    return {
        "job_id": job_id,
        "plots_generated": job.get("plots_generated", False),
        "plots_count": job.get("plots_count", 0),
        "plots_dir": job.get("plots_dir")
    }


@app.get("/plots/{job_id}/download")
async def download_plots(job_id: str):
    """
    Download all plots as a ZIP file
    """
    if job_id not in jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = jobs[job_id]
    
    if not job.get("plots_generated", False):
        raise HTTPException(status_code=400, detail="Plots not generated yet")
    
    plots_dir = Path(job.get("plots_dir"))
    if not plots_dir.exists():
        raise HTTPException(status_code=404, detail="Plots directory not found")
    
    # Create ZIP file
    import zipfile
    zip_path = RESULTS_DIR / f"{job_id}_plots.zip"
    
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for plot_file in plots_dir.glob("*.png"):
            zipf.write(plot_file, arcname=plot_file.name)
    
    return FileResponse(
        path=zip_path,
        filename=f"calibration_plots_{job_id}.zip",
        media_type="application/zip"
    )


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
