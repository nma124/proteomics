import { useState, useEffect } from 'react'
import axios from 'axios'
import Analytics from './Analytics'
import './App.css'

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000'

function App() {
  const [msFile, setMsFile] = useState(null)
  const [concFile, setConcFile] = useState(null)
  const [msFileId, setMsFileId] = useState(null)
  const [concFileId, setConcFileId] = useState(null)
  const [uploadingMs, setUploadingMs] = useState(false)
  const [uploadingConc, setUploadingConc] = useState(false)
  const [jobId, setJobId] = useState(null)
  const [jobStatus, setJobStatus] = useState(null)
  const [error, setError] = useState(null)
  const [polling, setPolling] = useState(false)

  // Poll job status
  useEffect(() => {
    if (!jobId || !polling) return

    const interval = setInterval(async () => {
      try {
        const response = await axios.get(`${API_BASE}/jobs/${jobId}`)
        setJobStatus(response.data)
        
        if (response.data.status === 'completed' || response.data.status === 'failed') {
          setPolling(false)
        }
      } catch (err) {
        console.error('Error polling job status:', err)
        setPolling(false)
      }
    }, 2000)

    return () => clearInterval(interval)
  }, [jobId, polling])

  const handleMsFileChange = (e) => {
    setMsFile(e.target.files[0])
    setMsFileId(null)
  }

  const handleConcFileChange = (e) => {
    setConcFile(e.target.files[0])
    setConcFileId(null)
  }

  const uploadMsFile = async () => {
    if (!msFile) return

    setUploadingMs(true)
    setError(null)
    
    const formData = new FormData()
    formData.append('file', msFile)

    try {
      const response = await axios.post(`${API_BASE}/upload/ms`, formData)
      setMsFileId(response.data.file_id)
    } catch (err) {
      setError(`Failed to upload MS file: ${err.response?.data?.detail || err.message}`)
    } finally {
      setUploadingMs(false)
    }
  }

  const uploadConcFile = async () => {
    if (!concFile) return

    setUploadingConc(true)
    setError(null)
    
    const formData = new FormData()
    formData.append('file', concFile)

    try {
      const response = await axios.post(`${API_BASE}/upload/concentration`, formData)
      setConcFileId(response.data.file_id)
    } catch (err) {
      setError(`Failed to upload concentration file: ${err.response?.data?.detail || err.message}`)
    } finally {
      setUploadingConc(false)
    }
  }

  const processData = async () => {
    if (!msFileId || !concFileId) return

    setError(null)
    
    try {
      const response = await axios.post(`${API_BASE}/process`, {
        ms_file_id: msFileId,
        concentration_file_id: concFileId,
        output_filename: 'results.csv'
      })
      
      setJobId(response.data.job_id)
      setJobStatus(response.data)
      setPolling(true)
    } catch (err) {
      setError(`Failed to start processing: ${err.response?.data?.detail || err.message}`)
    }
  }

  const downloadResults = () => {
    window.open(`${API_BASE}/download/${jobId}`, '_blank')
  }

  const reset = () => {
    setMsFile(null)
    setConcFile(null)
    setMsFileId(null)
    setConcFileId(null)
    setJobId(null)
    setJobStatus(null)
    setError(null)
    setPolling(false)
  }

  return (
    <div className="app">
      <header className="header">
        <h1>🧬 Proteomics PRM Data Processor</h1>
        <p>Unified pipeline for Skyline PRM data with automatic format detection</p>
      </header>

      <main className="main">
        {error && (
          <div className="error-box">
            <strong>❌ Error:</strong> {error}
          </div>
        )}

        {!jobId ? (
          <>
            {/* File Upload Section */}
            <section className="upload-section">
              <div className="upload-card">
                <h2>📁 Step 1: Upload MS Data</h2>
                <p className="help-text">Skyline PRM export CSV file</p>
                <input 
                  type="file" 
                  accept=".csv" 
                  onChange={handleMsFileChange}
                  disabled={uploadingMs}
                />
                {msFile && !msFileId && (
                  <button 
                    onClick={uploadMsFile} 
                    disabled={uploadingMs}
                    className="btn btn-primary"
                  >
                    {uploadingMs ? '⏳ Uploading...' : '⬆️ Upload MS File'}
                  </button>
                )}
                {msFileId && (
                  <div className="success-msg">✅ Uploaded: {msFile.name}</div>
                )}
              </div>

              <div className="upload-card">
                <h2>📊 Step 2: Upload Concentration Data</h2>
                <p className="help-text">Peptide dilution/concentration CSV file</p>
                <input 
                  type="file" 
                  accept=".csv" 
                  onChange={handleConcFileChange}
                  disabled={uploadingConc}
                />
                {concFile && !concFileId && (
                  <button 
                    onClick={uploadConcFile} 
                    disabled={uploadingConc}
                    className="btn btn-primary"
                  >
                    {uploadingConc ? '⏳ Uploading...' : '⬆️ Upload Concentration File'}
                  </button>
                )}
                {concFileId && (
                  <div className="success-msg">✅ Uploaded: {concFile.name}</div>
                )}
              </div>
            </section>

            {/* Process Button */}
            {msFileId && concFileId && (
              <section className="process-section">
                <button 
                  onClick={processData} 
                  className="btn btn-large btn-success"
                >
                  🚀 Process Data
                </button>
              </section>
            )}
          </>
        ) : (
          <>
            {/* Job Status Section */}
            <section className="status-section">
              <div className="status-card">
                <h2>📊 Processing Status</h2>
                
                <div className="status-info">
                  <div className="status-row">
                    <span className="label">Job ID:</span>
                    <code>{jobId}</code>
                  </div>
                  
                  <div className="status-row">
                    <span className="label">Status:</span>
                    <span className={`status-badge status-${jobStatus?.status}`}>
                      {jobStatus?.status?.toUpperCase()}
                    </span>
                  </div>

                  {jobStatus?.message && (
                    <div className="status-row">
                      <span className="label">Message:</span>
                      <span>{jobStatus.message}</span>
                    </div>
                  )}

                  {jobStatus?.progress !== undefined && (
                    <div className="progress-section">
                      <div className="progress-bar">
                        <div 
                          className="progress-fill" 
                          style={{ width: `${jobStatus.progress}%` }}
                        />
                      </div>
                      <span className="progress-text">{jobStatus.progress}%</span>
                    </div>
                  )}

                  {jobStatus?.result_stats && (
                    <div className="results-stats">
                      <h3>📈 Results Summary</h3>
                      <div className="stats-grid">
                        <div className="stat">
                          <div className="stat-value">{jobStatus.result_stats.rows}</div>
                          <div className="stat-label">Rows</div>
                        </div>
                        <div className="stat">
                          <div className="stat-value">{jobStatus.result_stats.columns}</div>
                          <div className="stat-label">Columns</div>
                        </div>
                        {jobStatus.result_stats.peptides && (
                          <div className="stat">
                            <div className="stat-value">{jobStatus.result_stats.peptides}</div>
                            <div className="stat-label">Peptides</div>
                          </div>
                        )}
                      </div>
                    </div>
                  )}
                </div>

                <div className="action-buttons">
                  {jobStatus?.status === 'completed' && (
                    <button 
                      onClick={downloadResults} 
                      className="btn btn-large btn-success"
                    >
                      💾 Download Results
                    </button>
                  )}
                  
                  <button 
                    onClick={reset} 
                    className="btn btn-secondary"
                  >
                    🔄 Process Another File
                  </button>
                </div>
              </div>
            </section>

            {/* Analytics Section */}
            {jobStatus?.status === 'completed' && (
              <section className="analytics-wrapper">
                <Analytics jobId={jobId} />
              </section>
            )}
          </>
        )}
      </main>

      <footer className="footer">
        <p>Automatic format detection • Heavy/Light paired • JPT single peptide</p>
      </footer>
    </div>
  )
}

export default App
