import { useState, useEffect } from 'react'
import axios from 'axios'
import './Analytics.css'

const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000'

function Analytics({ jobId }) {
  const [analytics, setAnalytics] = useState(null)
  const [loading, setLoading] = useState(true)
  const [error, setError] = useState(null)
  const [selectedPlot, setSelectedPlot] = useState(null)
  const [plotData, setPlotData] = useState(null)
  const [loadingPlot, setLoadingPlot] = useState(false)
  const [activeTab, setActiveTab] = useState('summary')
  const [generatingAll, setGeneratingAll] = useState(false)
  const [plotsStatus, setPlotsStatus] = useState({ plots_generated: false, plots_count: 0 })

  useEffect(() => {
    loadAnalytics()
    checkPlotsStatus()
  }, [jobId])

  const loadAnalytics = async () => {
    setLoading(true)
    setError(null)
    
    try {
      const response = await axios.get(`${API_BASE}/analytics/${jobId}`)
      setAnalytics(response.data)
      
      // Auto-select first plot category if available
      if (response.data.plot_categories && response.data.plot_categories.length > 0) {
        loadPlot(response.data.plot_categories[0])
      }
    } catch (err) {
      setError(`Failed to load analytics: ${err.response?.data?.detail || err.message}`)
    } finally {
      setLoading(false)
    }
  }

  const loadPlot = async (plotCategory) => {
    setSelectedPlot(plotCategory)
    setLoadingPlot(true)
    
    try {
      const response = await axios.get(`${API_BASE}/plots/${jobId}/${encodeURIComponent(plotCategory)}`)
      setPlotData(response.data)
    } catch (err) {
      console.error('Failed to load plot:', err)
      setPlotData(null)
    } finally {
      setLoadingPlot(false)
    }
  }

  const checkPlotsStatus = async () => {
    try {
      const response = await axios.get(`${API_BASE}/plots/${jobId}/status`)
      setPlotsStatus(response.data)
    } catch (err) {
      console.error('Failed to check plots status:', err)
    }
  }

  const generateAllPlots = async () => {
    setGeneratingAll(true)
    try {
      await axios.post(`${API_BASE}/plots/${jobId}/generate-all`)
      
      // Poll for completion
      const pollInterval = setInterval(async () => {
        const response = await axios.get(`${API_BASE}/plots/${jobId}/status`)
        setPlotsStatus(response.data)
        
        if (response.data.plots_generated) {
          clearInterval(pollInterval)
          setGeneratingAll(false)
        }
      }, 2000)
      
      // Stop polling after 5 minutes
      setTimeout(() => {
        clearInterval(pollInterval)
        setGeneratingAll(false)
      }, 300000)
      
    } catch (err) {
      console.error('Failed to generate plots:', err)
      setGeneratingAll(false)
    }
  }

  const downloadAllPlots = () => {
    window.open(`${API_BASE}/plots/${jobId}/download`, '_blank')
  }

  if (loading) {
    return <div className="analytics-loading">📊 Loading analytics...</div>
  }

  if (error) {
    return <div className="analytics-error">❌ {error}</div>
  }

  if (!analytics) {
    return <div className="analytics-error">No analytics data available</div>
  }

  return (
    <div className="analytics-container">
      <div className="analytics-tabs">
        <button 
          className={`tab ${activeTab === 'summary' ? 'active' : ''}`}
          onClick={() => setActiveTab('summary')}
        >
          📊 Summary
        </button>
        <button 
          className={`tab ${activeTab === 'quality' ? 'active' : ''}`}
          onClick={() => setActiveTab('quality')}
        >
          ✓ Quality Metrics
        </button>
        <button 
          className={`tab ${activeTab === 'peptides' ? 'active' : ''}`}
          onClick={() => setActiveTab('peptides')}
        >
          🧬 Per-Peptide Stats
        </button>
        <button 
          className={`tab ${activeTab === 'plots' ? 'active' : ''}`}
          onClick={() => setActiveTab('plots')}
        >
          📈 Calibration Plots
        </button>
      </div>

      <div className="analytics-content">
        {/* Summary Tab */}
        {activeTab === 'summary' && (
          <div className="summary-section">
            <h2>📊 Data Summary</h2>
            <div className="summary-grid">
              <div className="summary-card">
                <div className="summary-value">{analytics.summary.total_rows}</div>
                <div className="summary-label">Total Rows</div>
              </div>
              <div className="summary-card">
                <div className="summary-value">{analytics.summary.total_columns}</div>
                <div className="summary-label">Total Columns</div>
              </div>
              <div className="summary-card">
                <div className="summary-value">{analytics.summary.unique_peptides}</div>
                <div className="summary-label">Unique Peptides</div>
              </div>
              <div className="summary-card">
                <div className="summary-value">{analytics.summary.unique_fragments}</div>
                <div className="summary-label">Fragment Ions</div>
              </div>
              <div className="summary-card">
                <div className="summary-value">{analytics.summary.unique_plot_categories}</div>
                <div className="summary-label">Plot Categories</div>
              </div>
            </div>

            <h3>Available Data Columns</h3>
            <div className="columns-info">
              {analytics.available_columns.regression.length > 0 && (
                <div className="column-group">
                  <strong>Regression:</strong> {analytics.available_columns.regression.join(', ')}
                </div>
              )}
              {analytics.available_columns.aggregates.length > 0 && (
                <div className="column-group">
                  <strong>Aggregates:</strong> {analytics.available_columns.aggregates.join(', ')}
                </div>
              )}
              {analytics.available_columns.qtest.length > 0 && (
                <div className="column-group">
                  <strong>Q-test:</strong> {analytics.available_columns.qtest.join(', ')}
                </div>
              )}
            </div>
          </div>
        )}

        {/* Quality Metrics Tab */}
        {activeTab === 'quality' && (
          <div className="quality-section">
            <h2>✓ Regression Quality Metrics</h2>
            {analytics.quality_metrics && Object.keys(analytics.quality_metrics).length > 0 ? (
              <div className="quality-grid">
                <div className="quality-card">
                  <div className="quality-label">Mean R²</div>
                  <div className="quality-value">
                    {analytics.quality_metrics.r2_mean?.toFixed(4) || 'N/A'}
                  </div>
                </div>
                <div className="quality-card">
                  <div className="quality-label">Median R²</div>
                  <div className="quality-value">
                    {analytics.quality_metrics.r2_median?.toFixed(4) || 'N/A'}
                  </div>
                </div>
                <div className="quality-card">
                  <div className="quality-label">Min R²</div>
                  <div className="quality-value">
                    {analytics.quality_metrics.r2_min?.toFixed(4) || 'N/A'}
                  </div>
                </div>
                <div className="quality-card">
                  <div className="quality-label">Max R²</div>
                  <div className="quality-value">
                    {analytics.quality_metrics.r2_max?.toFixed(4) || 'N/A'}
                  </div>
                </div>
                <div className="quality-card good">
                  <div className="quality-label">Good Fits (R² &gt; 0.95)</div>
                  <div className="quality-value">{analytics.quality_metrics.r2_good_fits}</div>
                </div>
                <div className="quality-card poor">
                  <div className="quality-label">Poor Fits (R² &lt; 0.8)</div>
                  <div className="quality-value">{analytics.quality_metrics.r2_poor_fits}</div>
                </div>
              </div>
            ) : (
              <p>No quality metrics available</p>
            )}
          </div>
        )}

        {/* Per-Peptide Stats Tab */}
        {activeTab === 'peptides' && (
          <div className="peptides-section">
            <h2>🧬 Per-Peptide Statistics</h2>
            {analytics.peptide_stats && analytics.peptide_stats.length > 0 ? (
              <div className="table-container">
                <table className="peptide-table">
                  <thead>
                    <tr>
                      <th>Peptide</th>
                      <th>Measurements</th>
                      <th>Fragments</th>
                      <th>Mean R²</th>
                      <th>Std R²</th>
                      <th>Mean Gradient</th>
                    </tr>
                  </thead>
                  <tbody>
                    {analytics.peptide_stats.map((stat, idx) => (
                      <tr key={idx}>
                        <td className="peptide-name">{stat.peptide}</td>
                        <td>{stat.n_measurements}</td>
                        <td>{stat.n_fragments}</td>
                        <td className={stat.r2_mean > 0.95 ? 'good-r2' : stat.r2_mean < 0.8 ? 'poor-r2' : ''}>
                          {stat.r2_mean?.toFixed(4) || 'N/A'}
                        </td>
                        <td>{stat.r2_std?.toFixed(4) || 'N/A'}</td>
                        <td>{stat.gradient_mean?.toFixed(3) || 'N/A'}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <p>No peptide statistics available</p>
            )}
          </div>
        )}

        {/* Calibration Plots Tab */}
        {activeTab === 'plots' && (
          <div className="plots-section">
            <div className="plots-header">
              <h2>📈 Calibration Plots</h2>
              <div className="plots-actions">
                {!plotsStatus.plots_generated ? (
                  <button 
                    onClick={generateAllPlots}
                    disabled={generatingAll}
                    className="btn btn-generate"
                  >
                    {generatingAll ? '⏳ Generating All Plots...' : '🖼️ Generate All Plots'}
                  </button>
                ) : (
                  <div className="plots-ready">
                    <span className="plots-count">✅ {plotsStatus.plots_count} plots ready</span>
                    <button onClick={downloadAllPlots} className="btn btn-download">
                      📥 Download All Plots (ZIP)
                    </button>
                  </div>
                )}
              </div>
            </div>
            
            {analytics.plot_categories && analytics.plot_categories.length > 0 ? (
              <div className="plots-layout">
                <div className="plot-selector">
                  <h3>Select Plot Category:</h3>
                  <div className="plot-list">
                    {analytics.plot_categories.map((cat, idx) => (
                      <button
                        key={idx}
                        className={`plot-item ${selectedPlot === cat ? 'active' : ''}`}
                        onClick={() => loadPlot(cat)}
                      >
                        {cat}
                      </button>
                    ))}
                  </div>
                </div>

                <div className="plot-display">
                  {loadingPlot ? (
                    <div className="plot-loading">⏳ Loading plot...</div>
                  ) : plotData ? (
                    <div className="plot-result">
                      <h3>{plotData.plot_category}</h3>
                      <img 
                        src={plotData.image} 
                        alt={`Calibration plot for ${plotData.plot_category}`}
                        className="calibration-plot"
                      />
                      <div className="plot-stats">
                        <div className="stat-item">
                          <strong>R²:</strong> {plotData.stats.r2.toFixed(4)}
                        </div>
                        <div className="stat-item">
                          <strong>Slope:</strong> {plotData.stats.slope.toFixed(3)}
                        </div>
                        <div className="stat-item">
                          <strong>Intercept:</strong> {plotData.stats.intercept.toFixed(3)}
                        </div>
                        <div className="stat-item">
                          <strong>Points:</strong> {plotData.stats.n_points}
                        </div>
                      </div>
                    </div>
                  ) : (
                    <div className="plot-placeholder">
                      Select a plot category to view the calibration curve
                    </div>
                  )}
                </div>
              </div>
            ) : (
              <p>No plot categories available</p>
            )}
          </div>
        )}
      </div>
    </div>
  )
}

export default Analytics
