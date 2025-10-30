# Proteomics Web Interface Guide

A modern web interface for processing Skyline PRM data with automatic format detection.

## 🏗️ Architecture

```
proteomics/
├── api/                    # FastAPI Backend
│   ├── app.py             # Main API application
│   ├── uploads/           # Temporary file storage (auto-created)
│   └── results/           # Processed results (auto-created)
├── web/                    # React Frontend
│   ├── src/
│   │   ├── App.jsx        # Main React component
│   │   ├── App.css        # Styles
│   │   └── main.jsx       # Entry point
│   ├── index.html         # HTML template
│   ├── package.json       # Node dependencies
│   └── vite.config.js     # Vite configuration
└── main.py                # CLI (still works independently)
```

## 🚀 Quick Start

### 1. Install Backend Dependencies

```bash
# Install Python dependencies
pip install -r requirements.txt
```

### 2. Install Frontend Dependencies

```bash
# Navigate to web directory
cd web

# Install Node.js dependencies
npm install
```

### 3. Start the Backend Server

```bash
# From project root
python api/app.py

# Or with uvicorn directly
uvicorn api.app:app --reload --host 0.0.0.0 --port 8000
```

The API will be available at: http://localhost:8000
- API docs (Swagger): http://localhost:8000/docs
- Alternative docs (ReDoc): http://localhost:8000/redoc

### 4. Start the Frontend Server

```bash
# In a new terminal, from web/ directory
cd web
npm run dev
```

The web interface will be available at: http://localhost:3000

## 📋 Using the Web Interface

### Step-by-Step Workflow

1. **Upload MS Data File**
   - Click "Choose File" under Step 1
   - Select your Skyline PRM export CSV
   - Click "⬆️ Upload MS File"
   - Wait for ✅ confirmation

2. **Upload Concentration File**
   - Click "Choose File" under Step 2
   - Select your peptide concentration/dilution CSV
   - Click "⬆️ Upload Concentration File"
   - Wait for ✅ confirmation

3. **Process Data**
   - Once both files are uploaded, click "🚀 Process Data"
   - Monitor real-time progress with status updates
   - View processing statistics when complete

4. **Download Results**
   - Click "💾 Download Results" to get your processed CSV
   - Click "🔄 Process Another File" to start over

### Automatic Format Detection

The system automatically detects and processes:
- ✅ Heavy/Light paired peptides (area ratio analysis)
- ✅ JPT format (single peptide with conditions)
- ✅ Basic single peptide intensity data

No manual configuration needed!

## 🔧 API Endpoints

### File Upload
- `POST /upload/ms` - Upload mass spectrometry data
- `POST /upload/concentration` - Upload concentration data

### Processing
- `POST /process` - Start data processing job
- `GET /jobs/{job_id}` - Check job status
- `GET /jobs` - List all jobs

### Results
- `GET /download/{job_id}` - Download processed results
- `DELETE /jobs/{job_id}` - Delete job and results

## 🛠️ Development

### Backend Development

```bash
# Run with auto-reload
uvicorn api.app:app --reload

# Run with custom port
uvicorn api.app:app --port 8080
```

### Frontend Development

```bash
cd web

# Development server (hot reload)
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview
```

### CORS Configuration

The backend allows requests from:
- `http://localhost:3000` (Vite default)
- `http://localhost:5173` (Vite alternative)

To add more origins, edit `api/app.py`:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://your-domain.com"],
    ...
)
```

## 📊 Input File Requirements

### MS Data File (Skyline Export)
Required columns:
- `Peptide` - Peptide sequence
- `Protein` - Protein identifier
- `Replicate` - Sample replicate name
- `Precursor Mz` - Precursor mass-to-charge
- `Precursor Charge` - Precursor charge state
- `Product Mz` - Fragment ion m/z
- `Product Charge` - Fragment ion charge
- `Fragment Ion` - Fragment ion identifier
- `Area` - Peak area

### Concentration File
Required format:
- `Peptides` column with peptide sequences
- `D0 (ng/mL)`, `D1 (ng/mL)`, ..., `D7 (ng/mL)` columns with concentrations

## 🐛 Troubleshooting

### Backend Issues

**Port already in use:**
```bash
# Use different port
uvicorn api.app:app --port 8001
```

**Import errors:**
```bash
# Reinstall dependencies
pip install -r requirements.txt --force-reinstall
```

**File upload fails:**
- Check `api/uploads/` directory exists and is writable
- Verify file is valid CSV format
- Check file size (large files may timeout)

### Frontend Issues

**npm install fails:**
```bash
# Clear cache and reinstall
rm -rf node_modules package-lock.json
npm install
```

**Cannot connect to backend:**
- Ensure backend is running on port 8000
- Check CORS settings in `api/app.py`
- Verify proxy configuration in `web/vite.config.js`

**Build errors:**
```bash
# Update dependencies
npm update
```

## 🔒 Production Considerations

### Backend
1. **Use production ASGI server:**
   ```bash
   gunicorn api.app:app -w 4 -k uvicorn.workers.UvicornWorker
   ```

2. **Add authentication:**
   - Implement JWT tokens
   - Add API key validation
   - Use OAuth2

3. **Persistent storage:**
   - Replace in-memory `jobs` dict with Redis/PostgreSQL
   - Store files in S3/cloud storage
   - Add database for job history

4. **Security:**
   - Add rate limiting
   - Validate file contents
   - Scan for malicious uploads
   - Use HTTPS

### Frontend
1. **Build for production:**
   ```bash
   cd web
   npm run build
   ```

2. **Serve with nginx/Apache:**
   ```nginx
   server {
       listen 80;
       server_name your-domain.com;
       
       root /path/to/proteomics/web/dist;
       index index.html;
       
       location / {
           try_files $uri $uri/ /index.html;
       }
       
       location /api {
           proxy_pass http://localhost:8000;
       }
   }
   ```

3. **Environment variables:**
   - Create `.env` files for different environments
   - Configure API endpoints dynamically

## 📝 Example Usage

### Using cURL

```bash
# Upload MS file
curl -X POST http://localhost:8000/upload/ms \
  -F "file=@data/input/ms_data.csv"

# Response: {"file_id": "abc-123", ...}

# Upload concentration file
curl -X POST http://localhost:8000/upload/concentration \
  -F "file=@data/input/concentrations.csv"

# Response: {"file_id": "def-456", ...}

# Start processing
curl -X POST http://localhost:8000/process \
  -H "Content-Type: application/json" \
  -d '{"ms_file_id": "abc-123", "concentration_file_id": "def-456"}'

# Response: {"job_id": "xyz-789", "status": "pending", ...}

# Check status
curl http://localhost:8000/jobs/xyz-789

# Download results
curl -O http://localhost:8000/download/xyz-789
```

## 🎯 CLI vs Web Interface

Both interfaces use the same processing pipeline!

**CLI (Command Line):**
```bash
python main.py ms_data.csv concentrations.csv -o results.csv
```

**Web Interface:**
- Upload files through browser
- Monitor progress in real-time
- Download results with one click

Choose based on your workflow:
- CLI: Automation, scripting, batch processing
- Web: Interactive use, collaboration, visualization

## 🤝 Contributing

To add new features:

1. **Backend:** Add endpoints to `api/app.py`
2. **Frontend:** Update `web/src/App.jsx`
3. **Processing:** Modify `scripts/unified_processor.py`

## 📚 Additional Resources

- [FastAPI Documentation](https://fastapi.tiangolo.com/)
- [React Documentation](https://react.dev/)
- [Vite Documentation](https://vitejs.dev/)
- [Uvicorn Documentation](https://www.uvicorn.org/)

---

**Need help?** Check the main [README.md](README.md) or open an issue.
