# Proteomics Web Interface

A modern, user-friendly web interface for processing Skyline PRM (Parallel Reaction Monitoring) data with automatic format detection.

![Version](https://img.shields.io/badge/version-1.0.0-blue)
![React](https://img.shields.io/badge/react-18.2.0-61dafb)
![Vite](https://img.shields.io/badge/vite-5.0.0-646cff)

## ✨ Features

- 🎯 **Drag-and-drop file upload** - Upload MS data and concentration files
- 🔄 **Real-time progress tracking** - Monitor processing status with live updates
- 📊 **Automatic format detection** - Handles Heavy/Light paired, JPT, and single peptide formats
- 💾 **One-click download** - Get your processed results instantly
- 🎨 **Modern UI** - Beautiful, responsive design that works on all devices
- ⚡ **Fast processing** - Background job processing with status polling

## 🚀 Quick Start

### Prerequisites

- Node.js 16+ and npm
- Backend API running on `http://localhost:8000` (see `../api/README.md`)

### Installation

```bash
# Install dependencies
npm install
```

### Development

```bash
# Start development server with hot reload
npm run dev
```

The app will be available at **http://localhost:3000**

### Production Build

```bash
# Build for production
npm run build

# Preview production build
npm run preview
```

## 📖 Usage

### Step-by-Step Workflow

1. **Upload MS Data File**
   - Click "Choose File" under Step 1
   - Select your Skyline PRM export CSV file
   - Click "⬆️ Upload MS File"
   - Wait for the ✅ success confirmation

2. **Upload Concentration Data**
   - Click "Choose File" under Step 2
   - Select your peptide concentration/dilution CSV file
   - Click "⬆️ Upload Concentration File"
   - Wait for the ✅ success confirmation

3. **Process Data**
   - Once both files are uploaded, the "🚀 Process Data" button appears
   - Click to start processing
   - Watch real-time progress updates

4. **Download Results**
   - When processing completes, view summary statistics
   - Click "💾 Download Results" to get your CSV file
   - Click "🔄 Process Another File" to start a new analysis

### Supported File Formats

The system automatically detects and processes:

- ✅ **Heavy/Light paired peptides** - Area ratio analysis
- ✅ **JPT format** - Single peptide with multiple conditions
- ✅ **Basic format** - Single peptide intensity data

No manual configuration needed!

## 🏗️ Project Structure

```
web/
├── src/
│   ├── App.jsx          # Main application component
│   ├── App.css          # Styles and responsive design
│   └── main.jsx         # React entry point
├── index.html           # HTML template
├── package.json         # Dependencies and scripts
├── vite.config.js       # Vite configuration
└── README.md           # This file
```

## ⚙️ Configuration

### API Endpoint

By default, the app connects to `http://localhost:8000`. To change this, edit `src/App.jsx`:

```javascript
const API_BASE = 'http://your-api-server:8000'
```

### Proxy Configuration

The development server proxies `/api/*` requests to the backend. Configure in `vite.config.js`:

```javascript
server: {
  proxy: {
    '/api': {
      target: 'http://localhost:8000',
      changeOrigin: true
    }
  }
}
```

## 🎨 Customization

### Styling

All styles are in `src/App.css`. The design uses:
- CSS Grid for responsive layouts
- CSS custom properties for theming
- Smooth transitions and animations
- Mobile-first responsive design

### Color Scheme

Primary gradient:
```css
background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
```

To customize, search for these colors in `App.css`:
- Primary: `#667eea`
- Secondary: `#764ba2`
- Success: `#48bb78`
- Error: `#c33`

## 🔧 Development

### Available Scripts

```bash
# Start dev server
npm run dev

# Build for production
npm run build

# Preview production build
npm run preview
```

### Hot Module Replacement (HMR)

Vite provides instant HMR for a smooth development experience. Edit files in `src/` and see changes immediately.

### Debugging

Open browser DevTools to:
- View network requests to the API
- Check console logs for errors
- Inspect component state (use React DevTools)

## 🐛 Troubleshooting

### Cannot connect to backend

**Error:** `Network Error` when uploading files

**Solution:**
1. Ensure backend is running: `curl http://localhost:8000`
2. Check CORS settings in `../api/app.py`
3. Verify API_BASE URL in `src/App.jsx`

### npm install fails

**Solution:**
```bash
# Clear cache
rm -rf node_modules package-lock.json

# Reinstall
npm install
```

### Port 3000 already in use

**Solution:**
```bash
# Vite will automatically try port 3001, 3002, etc.
# Or specify a port:
npm run dev -- --port 3001
```

### Build errors

**Solution:**
```bash
# Update dependencies
npm update

# Clear Vite cache
rm -rf node_modules/.vite
npm run dev
```

## 📦 Dependencies

### Core
- **React 18.2** - UI library
- **Axios 1.6** - HTTP client for API requests

### Development
- **Vite 5.0** - Build tool and dev server
- **@vitejs/plugin-react** - React support for Vite

## 🚀 Deployment

### Static Hosting (Netlify, Vercel, GitHub Pages)

```bash
# Build
npm run build

# Deploy the dist/ directory
# Update API_BASE to point to production API
```

### Docker

```dockerfile
FROM node:18-alpine as build
WORKDIR /app
COPY package*.json ./
RUN npm ci
COPY . .
RUN npm run build

FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
```

### Environment Variables

For production, use environment variables:

```javascript
// src/config.js
export const API_BASE = import.meta.env.VITE_API_BASE || 'http://localhost:8000'
```

Then create `.env.production`:
```
VITE_API_BASE=https://api.yourdomain.com
```

## 🔒 Security Considerations

- Always use HTTPS in production
- Validate file types and sizes on both client and server
- Implement rate limiting for API requests
- Add authentication if handling sensitive data
- Sanitize user inputs

## 🤝 Contributing

To add features:

1. Create a new component in `src/components/`
2. Import and use in `App.jsx`
3. Add styles to `App.css`
4. Test thoroughly in both dev and prod builds

## 📝 Browser Support

- Chrome/Edge (latest)
- Firefox (latest)
- Safari (latest)
- Mobile browsers (iOS Safari, Chrome Mobile)

## 📚 Resources

- [React Documentation](https://react.dev/)
- [Vite Documentation](https://vitejs.dev/)
- [Axios Documentation](https://axios-http.com/)

## 📄 License

This project is part of the Proteomics PRM Processing toolkit. See main repository LICENSE.

---

**Need help?** Check the main [WEB_INTERFACE_GUIDE.md](../WEB_INTERFACE_GUIDE.md) or API documentation at [api/README.md](../api/README.md)
