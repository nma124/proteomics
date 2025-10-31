# Deployment Guide

Complete guide for deploying the Proteomics Web Interface to production.

## 🏗️ Architecture Overview

```
┌─────────────────┐         ┌──────────────────┐
│   Vercel        │         │   Railway/Render │
│   (Frontend)    │────────▶│   (Backend API)  │
│   React + Vite  │  HTTPS  │   FastAPI        │
└─────────────────┘         └──────────────────┘
```

## 🚀 Quick Deploy

### Step 1: Deploy Backend (Railway - Recommended)

**Why Railway?** Free tier, Python support, automatic deployments from Git.

1. **Sign up at [railway.app](https://railway.app)**

2. **Create New Project → Deploy from GitHub repo**

3. **Add environment variables:**
   ```
   PORT=8000
   PYTHONPATH=/app
   ```

4. **Create `Procfile` in project root:**
   ```bash
   echo "web: uvicorn api.app:app --host 0.0.0.0 --port \$PORT" > Procfile
   ```

5. **Create `railway.json`:**
   ```json
   {
     "build": {
       "builder": "NIXPACKS"
     },
     "deploy": {
       "startCommand": "uvicorn api.app:app --host 0.0.0.0 --port $PORT",
       "restartPolicyType": "ON_FAILURE",
       "restartPolicyMaxRetries": 10
     }
   }
   ```

6. **Deploy!** Railway will give you a URL like:
   `https://proteomics-api-production.up.railway.app`

### Step 2: Deploy Frontend (Vercel)

1. **Sign up at [vercel.com](https://vercel.com)**

2. **Import your Git repository**
   - Select the `web-interface` branch
   - Or merge to main first

3. **Configure Build Settings:**
   - **Framework Preset:** Vite
   - **Root Directory:** `web`
   - **Build Command:** `npm run build`
   - **Output Directory:** `dist`

4. **Add Environment Variable:**
   - Key: `VITE_API_BASE`
   - Value: Your Railway backend URL (e.g., `https://proteomics-api-production.up.railway.app`)

5. **Deploy!** Vercel will give you a URL like:
   `https://proteomics-web.vercel.app`

---

## 📋 Detailed Instructions

### Backend Deployment Options

#### Option 1: Railway (Recommended - Free Tier)

**Pros:**
- Free tier: 500 hours/month
- Automatic deployments from Git
- Easy environment variables
- Built-in monitoring

**Setup:**

1. Create `Procfile` in project root:
```bash
cd ~/Documents/projects/Bioinformatics/proteomics
cat > Procfile << 'EOF'
web: uvicorn api.app:app --host 0.0.0.0 --port $PORT
EOF
```

2. Create `railway.json`:
```json
{
  "build": {
    "builder": "NIXPACKS"
  },
  "deploy": {
    "startCommand": "uvicorn api.app:app --host 0.0.0.0 --port $PORT"
  }
}
```

3. Push to GitHub:
```bash
git add Procfile railway.json
git commit -m "Add Railway deployment config"
git push origin web-interface
```

4. On Railway dashboard:
   - New Project → Deploy from GitHub
   - Select your repo
   - Railway auto-detects Python and installs requirements.txt

5. Set environment variables in Railway dashboard:
   - `PORT=8000`
   - `PYTHONPATH=/app`

6. Copy your Railway URL (e.g., `https://proteomics-production.up.railway.app`)

---

#### Option 2: Render (Free Tier)

**Pros:**
- Free tier available
- Simple setup
- PostgreSQL database option

**Setup:**

1. Create `render.yaml` in project root:
```yaml
services:
  - type: web
    name: proteomics-api
    env: python
    buildCommand: pip install -r requirements.txt
    startCommand: uvicorn api.app:app --host 0.0.0.0 --port $PORT
    envVars:
      - key: PORT
        value: 8000
      - key: PYTHONPATH
        value: /opt/render/project/src
```

2. Go to [render.com](https://render.com)
3. New → Web Service → Connect GitHub repo
4. Render auto-detects settings
5. Deploy!

---

#### Option 3: Fly.io (Developer Friendly)

**Pros:**
- Good free tier
- Better for long-running processes
- Docker-based

**Setup:**

1. Install Fly CLI:
```bash
curl -L https://fly.io/install.sh | sh
```

2. Login and create app:
```bash
fly auth login
fly launch
```

3. Deploy:
```bash
fly deploy
```

---

### Frontend Deployment (Vercel)

#### Prerequisites

Make sure your code is pushed to GitHub:
```bash
git add .
git commit -m "Prepare for deployment"
git push origin web-interface
```

#### Vercel Setup

1. **Go to [vercel.com](https://vercel.com)** and sign in with GitHub

2. **Import Project:**
   - Click "Add New" → "Project"
   - Select your repository
   - Click "Import"

3. **Configure Project:**
   - **Framework Preset:** Vite
   - **Root Directory:** `web`
   - **Build Command:** `npm run build`
   - **Output Directory:** `dist`
   - **Install Command:** `npm install`

4. **Environment Variables:**
   Click "Environment Variables" and add:
   ```
   Name: VITE_API_BASE
   Value: https://your-backend-url.railway.app
   ```
   (Use your Railway/Render backend URL)

5. **Deploy:**
   - Click "Deploy"
   - Wait 2-3 minutes
   - You'll get a URL like `https://proteomics-web.vercel.app`

#### Custom Domain (Optional)

In Vercel dashboard:
- Go to Project Settings → Domains
- Add your custom domain
- Update DNS records as instructed

---

## 🔧 Configuration Files

### For Backend

Create these in project root:

**Procfile:**
```
web: uvicorn api.app:app --host 0.0.0.0 --port $PORT
```

**railway.json:**
```json
{
  "build": {
    "builder": "NIXPACKS"
  },
  "deploy": {
    "startCommand": "uvicorn api.app:app --host 0.0.0.0 --port $PORT",
    "restartPolicyType": "ON_FAILURE"
  }
}
```

### For Frontend

Already created: `web/vercel.json`

Create `.env.production` in `web/`:
```
VITE_API_BASE=https://your-backend-url.railway.app
```

---

## 🔒 Security Checklist

Before deploying to production:

### Backend
- [ ] Update CORS origins in `api/app.py`:
  ```python
  allow_origins=[
      "https://proteomics-web.vercel.app",  # Your Vercel URL
      "http://localhost:3000"  # Keep for local dev
  ]
  ```
- [ ] Add rate limiting
- [ ] Add authentication (JWT)
- [ ] Enable HTTPS only
- [ ] Set file size limits
- [ ] Add input validation
- [ ] Use production ASGI server (Gunicorn)

### Frontend
- [ ] Update API_BASE in `web/src/App.jsx` via env variable
- [ ] Remove console.log statements
- [ ] Enable production build optimizations
- [ ] Add error tracking (Sentry)

---

## 🧪 Testing Deployment

### Test Backend

```bash
# Check health
curl https://your-backend-url.railway.app/

# Check interactive docs
open https://your-backend-url.railway.app/docs
```

### Test Frontend

1. Open your Vercel URL
2. Upload test files
3. Check browser console for errors
4. Verify API calls work
5. Test complete workflow

---

## 🐛 Troubleshooting

### Backend Issues

**Problem:** 502 Bad Gateway
- Check Railway logs
- Verify PORT environment variable
- Ensure requirements.txt includes all dependencies

**Problem:** Import errors
- Add `PYTHONPATH=/app` environment variable
- Check file structure matches local

**Problem:** File upload fails
- Check file size limits
- Verify disk space
- Ensure `api/uploads/` directory exists

### Frontend Issues

**Problem:** API requests fail (CORS)
- Update CORS origins in backend
- Verify API_BASE URL is correct
- Check browser console for errors

**Problem:** Build fails
- Clear node_modules: `rm -rf node_modules package-lock.json`
- Reinstall: `npm install`
- Check build locally: `npm run build`

**Problem:** Environment variables not working
- Restart Vercel deployment
- Verify variable name: `VITE_API_BASE`
- Check it's set in Vercel dashboard

---

## 💰 Cost Estimates

### Free Tier (Good for Development/Testing)

| Service | Free Tier | Limits |
|---------|-----------|--------|
| **Railway** | 500 hours/month | $5/month credit |
| **Vercel** | Unlimited | 100GB bandwidth |
| **Render** | 750 hours/month | Sleeps after 15min inactivity |

### Paid Plans (For Production)

| Service | Cost | Features |
|---------|------|----------|
| **Railway** | ~$5-20/month | Always on, more resources |
| **Vercel** | Free-$20/month | Custom domains, analytics |
| **Render** | $7/month | Always on web service |

---

## 📊 Monitoring

### Backend Monitoring

Railway provides built-in monitoring:
- CPU/Memory usage
- Request logs
- Error tracking

### Frontend Monitoring

Vercel provides:
- Deployment logs
- Function invocations
- Web analytics (paid)

### Add External Monitoring (Optional)

**Sentry** (Error Tracking):
```bash
npm install @sentry/react
```

**Uptime Monitoring:**
- UptimeRobot (free)
- Pingdom
- Better Uptime

---

## 🔄 CI/CD (Automatic Deployments)

Both Vercel and Railway support automatic deployments:

1. **Push to GitHub**
   ```bash
   git push origin web-interface
   ```

2. **Auto-deploy triggers:**
   - Backend: Railway auto-deploys
   - Frontend: Vercel auto-deploys

3. **Preview deployments:**
   - Vercel creates preview URLs for PRs
   - Railway can deploy from branches

---

## 📚 Next Steps

After deployment:

1. **Test thoroughly** with real data
2. **Monitor logs** for errors
3. **Set up alerts** for downtime
4. **Add authentication** for production use
5. **Configure custom domain**
6. **Enable HTTPS** (automatic on Vercel/Railway)
7. **Add database** for persistent job storage
8. **Implement caching** for better performance

---

## 🆘 Support

- **Railway Docs:** https://docs.railway.app/
- **Vercel Docs:** https://vercel.com/docs
- **FastAPI Deployment:** https://fastapi.tiangolo.com/deployment/
- **Vite Deployment:** https://vitejs.dev/guide/static-deploy.html

---

**Ready to deploy?** Start with Railway for backend, then Vercel for frontend!
