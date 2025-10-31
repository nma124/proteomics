# Vercel Deployment Fix

## The Issue
Vercel couldn't find Vite modules because the project structure has the frontend in a `web/` subdirectory.

## ✅ Solution: Configure Vercel Settings

When deploying to Vercel, use these **exact settings**:

### Framework Preset
- **Framework:** `Vite`

### Build & Development Settings
- **Root Directory:** `web` ⚠️ **IMPORTANT**
- **Build Command:** `npm run build` (leave as default)
- **Output Directory:** `dist` (leave as default)
- **Install Command:** `npm install` (leave as default)

### Environment Variables
Add this variable:
- **Name:** `VITE_API_BASE`
- **Value:** Your Railway backend URL (e.g., `https://your-app.up.railway.app`)

## 🚀 Step-by-Step Deployment

1. **Go to Vercel Dashboard**
   - Visit [vercel.com](https://vercel.com)
   - Click "Add New" → "Project"

2. **Import Repository**
   - Select your GitHub repository
   - Click "Import"

3. **Configure Project** ⚠️ Critical Step
   ```
   Framework Preset: Vite
   Root Directory: web          ← MUST BE SET!
   Build Command: npm run build
   Output Directory: dist
   Install Command: npm install
   ```

4. **Add Environment Variable**
   - Click "Environment Variables"
   - Add:
     - Name: `VITE_API_BASE`
     - Value: `https://your-railway-backend.up.railway.app`

5. **Deploy**
   - Click "Deploy"
   - Wait 2-3 minutes
   - Done! 🎉

## 🐛 Still Having Issues?

### Option 1: Clear Build Cache
In Vercel dashboard:
- Go to Settings → General
- Scroll to "Build & Development Settings"
- Click "Clear Cache"
- Redeploy

### Option 2: Manual Override
If auto-detection fails, manually override:
```json
{
  "buildCommand": "npm run build",
  "outputDirectory": "dist", 
  "installCommand": "npm install",
  "framework": "vite"
}
```

### Option 3: Check Node Version
In Vercel Settings → General → Node.js Version:
- Set to: `18.x` or `20.x`

## ✅ Verification

After deployment:
1. Visit your Vercel URL
2. Open browser console (F12)
3. Check for errors
4. Try uploading a file
5. Verify API calls work

## 📝 Project Structure
```
proteomics/
├── web/              ← Set this as Root Directory in Vercel
│   ├── src/
│   ├── dist/         ← Build output
│   ├── package.json
│   └── vite.config.js
└── api/              ← Backend (not deployed to Vercel)
```

## 🔗 Helpful Links

- [Vercel Vite Guide](https://vercel.com/docs/frameworks/vite)
- [Vercel Environment Variables](https://vercel.com/docs/projects/environment-variables)
- [Vercel Troubleshooting](https://vercel.com/docs/errors)

---

**TL;DR:** Set **Root Directory** to `web` in Vercel project settings!
