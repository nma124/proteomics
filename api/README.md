# Proteomics API Documentation

FastAPI backend for processing Skyline PRM data with automatic format detection.

![Python](https://img.shields.io/badge/python-3.7+-blue)
![FastAPI](https://img.shields.io/badge/fastapi-0.104+-green)
![License](https://img.shields.io/badge/license-MIT-blue)

## 📋 Overview

RESTful API that provides:
- File upload endpoints for MS data and concentration files
- Background job processing with status tracking
- Automatic format detection (Heavy/Light pairs, JPT, single peptide)
- Results download and job management

## 🚀 Quick Start

### Installation

```bash
# Install dependencies
pip install fastapi "uvicorn[standard]" python-multipart

# Or from project root
pip install -r requirements.txt
```

### Run the Server

```bash
# Development mode (auto-reload)
uvicorn api.app:app --reload --host 0.0.0.0 --port 8000

# Production mode
python api/app.py
```

The API will be available at:
- **API:** http://localhost:8000
- **Interactive docs (Swagger):** http://localhost:8000/docs
- **Alternative docs (ReDoc):** http://localhost:8000/redoc

## 📡 API Endpoints

### Health Check

#### `GET /`
Check if the API is online.

**Response:**
```json
{
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
```

---

### File Upload

#### `POST /upload/ms`
Upload mass spectrometry data file (Skyline PRM export CSV).

**Request:**
- **Content-Type:** `multipart/form-data`
- **Body:** `file` - CSV file

**Response:**
```json
{
  "file_id": "abc-123-def-456",
  "filename": "ms_data.csv",
  "size": 1048576,
  "type": "mass_spectrometry"
}
```

**Example:**
```bash
curl -X POST http://localhost:8000/upload/ms \
  -F "file=@data/ms_data.csv"
```

---

#### `POST /upload/concentration`
Upload peptide concentration/dilution CSV file.

**Request:**
- **Content-Type:** `multipart/form-data`
- **Body:** `file` - CSV file

**Response:**
```json
{
  "file_id": "xyz-789-uvw-012",
  "filename": "concentrations.csv",
  "size": 2048,
  "type": "concentration"
}
```

**Example:**
```bash
curl -X POST http://localhost:8000/upload/concentration \
  -F "file=@data/concentrations.csv"
```

---

### Data Processing

#### `POST /process`
Start processing uploaded data files.

**Request:**
```json
{
  "ms_file_id": "abc-123-def-456",
  "concentration_file_id": "xyz-789-uvw-012",
  "output_filename": "results.csv"  // optional
}
```

**Response:**
```json
{
  "job_id": "job-111-222-333",
  "status": "pending",
  "message": "Job created, starting processing...",
  "progress": 0,
  "created_at": "2024-01-15T10:30:00.000Z",
  "result_file": null,
  "completed_at": null
}
```

**Example:**
```bash
curl -X POST http://localhost:8000/process \
  -H "Content-Type: application/json" \
  -d '{
    "ms_file_id": "abc-123-def-456",
    "concentration_file_id": "xyz-789-uvw-012"
  }'
```

---

### Job Management

#### `GET /jobs/{job_id}`
Get status and details of a processing job.

**Response:**
```json
{
  "job_id": "job-111-222-333",
  "status": "completed",  // "pending", "processing", "completed", "failed"
  "message": "Processing complete! Generated 220 rows × 45 columns",
  "progress": 100,
  "created_at": "2024-01-15T10:30:00.000Z",
  "completed_at": "2024-01-15T10:31:30.000Z",
  "result_file": "/path/to/results.csv",
  "result_stats": {
    "rows": 220,
    "columns": 45,
    "peptides": 1
  }
}
```

**Status Values:**
- `pending` - Job created, waiting to start
- `processing` - Currently processing data
- `completed` - Processing finished successfully
- `failed` - Processing failed with error

**Example:**
```bash
curl http://localhost:8000/jobs/job-111-222-333
```

---

#### `GET /jobs`
List all jobs in the system.

**Response:**
```json
{
  "jobs": [
    {
      "job_id": "job-111-222-333",
      "status": "completed",
      "created_at": "2024-01-15T10:30:00.000Z",
      ...
    },
    {
      "job_id": "job-444-555-666",
      "status": "processing",
      "created_at": "2024-01-15T10:35:00.000Z",
      ...
    }
  ]
}
```

**Example:**
```bash
curl http://localhost:8000/jobs
```

---

#### `DELETE /jobs/{job_id}`
Delete a job and its associated files.

**Response:**
```json
{
  "message": "Job deleted successfully"
}
```

**Example:**
```bash
curl -X DELETE http://localhost:8000/jobs/job-111-222-333
```

---

### Results Download

#### `GET /download/{job_id}`
Download processed results CSV file.

**Response:**
- **Content-Type:** `text/csv`
- **Body:** CSV file content

**Example:**
```bash
curl -O http://localhost:8000/download/job-111-222-333
```

---

## 🏗️ Architecture

```
api/
├── app.py              # Main FastAPI application
├── uploads/           # Temporary uploaded files (auto-created)
├── results/           # Processed result files (auto-created)
└── README.md         # This file
```

### Data Flow

1. **Upload Phase**
   - Files uploaded via `/upload/ms` and `/upload/concentration`
   - Stored in `api/uploads/` with unique IDs
   - File IDs returned to client

2. **Processing Phase**
   - Client sends file IDs to `/process`
   - Background task created
   - Unified processor invoked (automatic format detection)
   - Job status updated in real-time

3. **Results Phase**
   - Processed data saved to `api/results/`
   - Statistics computed and stored
   - Client polls `/jobs/{job_id}` for status
   - Download via `/download/{job_id}` when complete

### Job Storage

Currently uses in-memory dictionary for job tracking. For production, consider:
- **Redis** - Fast, persistent job queue
- **PostgreSQL** - Full relational database
- **MongoDB** - Document-based storage

## ⚙️ Configuration

### CORS Settings

The API allows requests from:
- `http://localhost:3000` (Vite default)
- `http://localhost:5173` (Vite alternative)

To add origins, edit `app.py`:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://your-domain.com",
        "https://your-production-domain.com"
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### Port Configuration

```bash
# Default port 8000
uvicorn api.app:app

# Custom port
uvicorn api.app:app --port 8080

# Custom host and port
uvicorn api.app:app --host 0.0.0.0 --port 8080
```

### File Storage

By default, files are stored in:
- `api/uploads/` - Uploaded files
- `api/results/` - Processed results

To change, edit `app.py`:

```python
UPLOAD_DIR = Path("custom/upload/path")
RESULTS_DIR = Path("custom/results/path")
```

## 🔒 Security

### Current Implementation
- File type validation (CSV only)
- Unique file IDs (UUID4)
- CORS protection

### Production Recommendations

1. **Authentication**
   ```python
   from fastapi.security import HTTPBearer
   
   security = HTTPBearer()
   
   @app.post("/process")
   async def process_data(
       request: ProcessRequest,
       credentials: HTTPAuthorizationCredentials = Depends(security)
   ):
       # Validate token
       ...
   ```

2. **Rate Limiting**
   ```python
   from slowapi import Limiter
   
   limiter = Limiter(key_func=get_remote_address)
   
   @app.post("/upload/ms")
   @limiter.limit("10/minute")
   async def upload_ms_file(...):
       ...
   ```

3. **File Size Limits**
   ```python
   from fastapi import File, UploadFile
   
   MAX_FILE_SIZE = 100 * 1024 * 1024  # 100MB
   
   @app.post("/upload/ms")
   async def upload_ms_file(file: UploadFile = File(..., max_size=MAX_FILE_SIZE)):
       ...
   ```

4. **Input Validation**
   - Validate CSV structure
   - Check for malicious content
   - Sanitize filenames

5. **HTTPS Only**
   ```bash
   uvicorn api.app:app --ssl-keyfile=key.pem --ssl-certfile=cert.pem
   ```

## 🧪 Testing

### Manual Testing with cURL

```bash
# 1. Upload MS file
MS_RESPONSE=$(curl -s -X POST http://localhost:8000/upload/ms \
  -F "file=@test_data.csv")
MS_ID=$(echo $MS_RESPONSE | jq -r '.file_id')

# 2. Upload concentration file
CONC_RESPONSE=$(curl -s -X POST http://localhost:8000/upload/concentration \
  -F "file=@test_conc.csv")
CONC_ID=$(echo $CONC_RESPONSE | jq -r '.file_id')

# 3. Start processing
JOB_RESPONSE=$(curl -s -X POST http://localhost:8000/process \
  -H "Content-Type: application/json" \
  -d "{\"ms_file_id\": \"$MS_ID\", \"concentration_file_id\": \"$CONC_ID\"}")
JOB_ID=$(echo $JOB_RESPONSE | jq -r '.job_id')

# 4. Poll status
while true; do
  STATUS=$(curl -s http://localhost:8000/jobs/$JOB_ID | jq -r '.status')
  echo "Status: $STATUS"
  [ "$STATUS" = "completed" ] || [ "$STATUS" = "failed" ] && break
  sleep 2
done

# 5. Download results
curl -O http://localhost:8000/download/$JOB_ID
```

### Automated Testing (pytest)

```python
# tests/test_api.py
import pytest
from fastapi.testclient import TestClient
from api.app import app

client = TestClient(app)

def test_health_check():
    response = client.get("/")
    assert response.status_code == 200
    assert response.json()["status"] == "online"

def test_upload_ms_file():
    with open("test_data.csv", "rb") as f:
        response = client.post(
            "/upload/ms",
            files={"file": ("test.csv", f, "text/csv")}
        )
    assert response.status_code == 200
    assert "file_id" in response.json()
```

Run tests:
```bash
pytest tests/test_api.py -v
```

## 🐛 Troubleshooting

### Import Errors

**Error:** `ModuleNotFoundError: No module named 'scripts'`

**Solution:**
```bash
# Run from project root
cd /path/to/proteomics
python api/app.py

# Or add to PYTHONPATH
export PYTHONPATH=/path/to/proteomics:$PYTHONPATH
```

### Port Already in Use

**Error:** `Address already in use`

**Solution:**
```bash
# Find process using port 8000
lsof -ti:8000

# Kill process
kill $(lsof -ti:8000)

# Or use different port
uvicorn api.app:app --port 8001
```

### File Upload Fails

**Error:** `Upload failed: file not found`

**Solution:**
1. Check `api/uploads/` directory exists and is writable
2. Verify file path in request
3. Check disk space
4. Verify CSV format

### Processing Fails

Check logs for detailed error messages. Common issues:
- Missing columns in CSV
- Incompatible data format
- Memory issues with large files

## 📊 Performance

### Optimization Tips

1. **Use async/await for I/O operations**
2. **Implement caching for repeated processing**
3. **Add database indexing for job lookup**
4. **Use CDN for static files**
5. **Enable gzip compression**

### Monitoring

```python
import time
from fastapi import Request

@app.middleware("http")
async def add_process_time_header(request: Request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(process_time)
    return response
```

## 🚀 Deployment

### Using Gunicorn (Production)

```bash
pip install gunicorn

# Run with 4 workers
gunicorn api.app:app \
  -w 4 \
  -k uvicorn.workers.UvicornWorker \
  --bind 0.0.0.0:8000
```

### Docker

```dockerfile
FROM python:3.9-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000"]
```

Build and run:
```bash
docker build -t proteomics-api .
docker run -p 8000:8000 proteomics-api
```

### Systemd Service

```ini
# /etc/systemd/system/proteomics-api.service
[Unit]
Description=Proteomics API
After=network.target

[Service]
User=www-data
Group=www-data
WorkingDirectory=/opt/proteomics
Environment="PATH=/opt/proteomics/venv/bin"
ExecStart=/opt/proteomics/venv/bin/uvicorn api.app:app --host 0.0.0.0 --port 8000

[Install]
WantedBy=multi-user.target
```

Enable and start:
```bash
sudo systemctl enable proteomics-api
sudo systemctl start proteomics-api
sudo systemctl status proteomics-api
```

## 📚 API Client Examples

### Python

```python
import requests

# Upload files
with open('ms_data.csv', 'rb') as f:
    ms_response = requests.post(
        'http://localhost:8000/upload/ms',
        files={'file': f}
    )
ms_id = ms_response.json()['file_id']

with open('concentrations.csv', 'rb') as f:
    conc_response = requests.post(
        'http://localhost:8000/upload/concentration',
        files={'file': f}
    )
conc_id = conc_response.json()['file_id']

# Process
process_response = requests.post(
    'http://localhost:8000/process',
    json={
        'ms_file_id': ms_id,
        'concentration_file_id': conc_id
    }
)
job_id = process_response.json()['job_id']

# Poll status
while True:
    status_response = requests.get(f'http://localhost:8000/jobs/{job_id}')
    status = status_response.json()['status']
    if status in ['completed', 'failed']:
        break
    time.sleep(2)

# Download
if status == 'completed':
    results = requests.get(f'http://localhost:8000/download/{job_id}')
    with open('results.csv', 'wb') as f:
        f.write(results.content)
```

### JavaScript

```javascript
// Upload and process
async function processProteomicsData(msFile, concFile) {
  // Upload MS file
  const msFormData = new FormData();
  msFormData.append('file', msFile);
  const msResponse = await fetch('http://localhost:8000/upload/ms', {
    method: 'POST',
    body: msFormData
  });
  const { file_id: msId } = await msResponse.json();

  // Upload concentration file
  const concFormData = new FormData();
  concFormData.append('file', concFile);
  const concResponse = await fetch('http://localhost:8000/upload/concentration', {
    method: 'POST',
    body: concFormData
  });
  const { file_id: concId } = await concResponse.json();

  // Start processing
  const processResponse = await fetch('http://localhost:8000/process', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      ms_file_id: msId,
      concentration_file_id: concId
    })
  });
  const { job_id } = await processResponse.json();

  // Poll status
  return new Promise((resolve, reject) => {
    const interval = setInterval(async () => {
      const statusResponse = await fetch(`http://localhost:8000/jobs/${job_id}`);
      const { status } = await statusResponse.json();
      
      if (status === 'completed') {
        clearInterval(interval);
        resolve(job_id);
      } else if (status === 'failed') {
        clearInterval(interval);
        reject(new Error('Processing failed'));
      }
    }, 2000);
  });
}
```

## 📖 Additional Resources

- [FastAPI Documentation](https://fastapi.tiangolo.com/)
- [Uvicorn Documentation](https://www.uvicorn.org/)
- [Pydantic Documentation](https://docs.pydantic.dev/)
- [OpenAPI Specification](https://swagger.io/specification/)

## 🤝 Contributing

To extend the API:

1. Add new endpoints to `app.py`
2. Update Pydantic models for request/response
3. Add error handling
4. Update this documentation
5. Add tests

## 📄 License

Part of the Proteomics PRM Processing toolkit. See main repository LICENSE.

---

**Need help?** Check the [Web Interface Guide](../WEB_INTERFACE_GUIDE.md) or [Frontend Documentation](../web/README.md)
