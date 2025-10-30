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
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import Optional, Dict
import sys
from pathlib import Path
import uuid
import shutil
import json
from datetime import datetime

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


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
