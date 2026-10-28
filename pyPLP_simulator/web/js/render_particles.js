/**
 * pyPLP Simulator - dNTP Simulation, Particles & Fluorometer Graph
 */

// Mini Real-Time Fluorometer Graph
function drawFluorometerGraph() {
  const fc = (typeof fCanvas !== 'undefined' && fCanvas) ? fCanvas : document.getElementById('fluoroCanvas');
  if (!fc) return;
  const fContext = (typeof fCtx !== 'undefined' && fCtx) ? fCtx : fc.getContext('2d');
  if (!fContext) return;

  fContext.clearRect(0, 0, fc.width, fc.height);
  const w = fc.width;
  const h = fc.height;

  fContext.strokeStyle = "rgba(234, 179, 8, 0.3)";
  fContext.lineWidth = 1;
  fContext.beginPath();
  fContext.moveTo(0, h * 0.5);
  fContext.lineTo(w, h * 0.5);
  fContext.stroke();

  if (state.fluoroHistory.length > 1) {
    fContext.strokeStyle = "#eab308";
    fContext.lineWidth = 2;
    fContext.shadowColor = "#eab308";
    fContext.shadowBlur = 8;
    fContext.beginPath();

    const step = w / Math.max(state.fluoroHistory.length - 1, 1);
    for (let i = 0; i < state.fluoroHistory.length; i++) {
      const val = state.fluoroHistory[i];
      const x = i * step;
      const y = h - (val / 100) * (h - 8) - 4;
      if (i === 0) fContext.moveTo(x, y);
      else fContext.lineTo(x, y);
    }
    fContext.stroke();
  }
}

// --- 8.6 DYNAMIC dNTP MOLECULAR POOL & RECRUITMENT SIMULATOR ---
function updateAndDrawDntps(targetCtx, w, h, dt, phi29Pos, isSynthesizing, neededDntp) {
  const c = targetCtx || ctx;
  if (!state.dntpPool) state.dntpPool = [];
  if (!state.releasedPpi) state.releasedPpi = [];

  const types = ['dATP', 'dTTP', 'dCTP', 'dGTP'];
  const baseMap = { 'dATP': 'A', 'dTTP': 'T', 'dCTP': 'C', 'dGTP': 'G' };

  // Freeze particle motion and recruitment when simulation is paused
  const effectiveDt = (state.isPlaying && dt > 0) ? dt : 0;

  // Maintain floating dNTP molecules in the chamber around Phi29
  while (state.dntpPool.length < 24) {
    let type;
    if (isSynthesizing && neededDntp && Math.random() < 0.45) {
      type = neededDntp; // Ensure regular supply of cognate dNTP in solution
    } else {
      type = types[Math.floor(Math.random() * types.length)];
    }
    const angle = Math.random() * Math.PI * 2;
    const dist = 120 + Math.random() * 220;
    state.dntpPool.push({
      x: phi29Pos.x + Math.cos(angle) * dist,
      y: phi29Pos.y + Math.sin(angle) * dist,
      vx: (Math.random() - 0.5) * 30,
      vy: (Math.random() - 0.5) * 30,
      type: type,
      base: baseMap[type]
    });
  }

  c.save();
  for (let i = state.dntpPool.length - 1; i >= 0; i--) {
    const d = state.dntpPool[i];
    const isCognate = (neededDntp ? d.type === neededDntp : true);

    if (effectiveDt > 0 && isSynthesizing && phi29Pos) {
      const dx = phi29Pos.x - d.x;
      const dy = phi29Pos.y - d.y;
      const dist = Math.hypot(dx, dy) || 1;

      if (isCognate) {
        // COGNATE NUCLEOTIDE: Specifically recruited into Phi29 active site!
        if (dist < 280) {
          const speed = 80 + (280 - dist) * 1.3;
          d.vx += (dx / dist) * speed * effectiveDt * 2.8;
          d.vy += (dy / dist) * speed * effectiveDt * 2.8;
          d.vx *= 0.91;
          d.vy *= 0.91;
        }

        if (dist < 22) {
          // Catalytic incorporation into growing 3'-OH end! Cleave triphosphate -> spawn PPi
          state.releasedPpi.push({
            x: phi29Pos.x + (Math.random() - 0.5) * 14,
            y: phi29Pos.y - 12 - Math.random() * 8,
            vx: (Math.random() - 0.5) * 35 + 20,
            vy: -30 - Math.random() * 25,
            alpha: 1.0
          });
          state.dntpPool.splice(i, 1);
          continue;
        }
      } else {
        // NON-COGNATE NUCLEOTIDE: Polymerase proofreading / steric cleft deflects mismatch away
        if (dist < 55) {
          const repel = 65 * (1 - dist / 55);
          d.vx -= (dx / dist) * repel * effectiveDt * 3.0;
          d.vy -= (dy / dist) * repel * effectiveDt * 3.0;
        }
      }
    }

    if (effectiveDt > 0) {
      d.x += d.vx * effectiveDt;
      d.y += d.vy * effectiveDt;

      // Keep dNTPs gently contained in solution near the active chamber
      const distFromPhi = Math.hypot(d.x - phi29Pos.x, d.y - phi29Pos.y);
      if (distFromPhi > 400) {
        const angle = Math.random() * Math.PI * 2;
        const dist = 140 + Math.random() * 200;
        d.x = phi29Pos.x + Math.cos(angle) * dist;
        d.y = phi29Pos.y + Math.sin(angle) * dist;
        d.vx = (Math.random() - 0.5) * 28;
        d.vy = (Math.random() - 0.5) * 28;
      }
    }

    // Draw dNTP capsule with base color
    const style = getBaseStyle(d.base);
    c.save();
    c.fillStyle = "rgba(15, 23, 42, 0.92)";
    c.strokeStyle = style.border;
    c.lineWidth = (isCognate && isSynthesizing) ? 2.0 : 1.3;
    c.shadowColor = style.border;
    c.shadowBlur = (isCognate && isSynthesizing) ? 10 : 5;

    const capW = 28;
    const capH = 13;
    c.beginPath();
    c.roundRect(d.x - capW * 0.5, d.y - capH * 0.5, capW, capH, 6);
    c.fill();
    c.stroke();

    c.fillStyle = style.border;
    c.font = "bold 8px 'JetBrains Mono', monospace";
    c.textAlign = "center";
    c.textBaseline = "middle";
    c.fillText(d.type, d.x, d.y);

    // Triphosphate group (amber phosphate beads)
    c.fillStyle = "#f59e0b";
    c.beginPath();
    c.arc(d.x + capW * 0.5 + 3, d.y, 1.8, 0, Math.PI * 2);
    c.arc(d.x + capW * 0.5 + 7, d.y, 1.8, 0, Math.PI * 2);
    c.arc(d.x + capW * 0.5 + 11, d.y, 1.8, 0, Math.PI * 2);
    c.fill();

    c.restore();
  }

  // Draw released PPi (inorganic pyrophosphate) drifting away
  for (let i = state.releasedPpi.length - 1; i >= 0; i--) {
    const ppi = state.releasedPpi[i];
    if (effectiveDt > 0) {
      ppi.x += ppi.vx * effectiveDt;
      ppi.y += ppi.vy * effectiveDt;
      ppi.alpha -= effectiveDt * 1.1;
    }

    if (ppi.alpha <= 0) {
      state.releasedPpi.splice(i, 1);
      continue;
    }

    c.save();
    c.globalAlpha = Math.max(0, ppi.alpha);
    c.fillStyle = "#fbbf24";
    c.shadowColor = "#f59e0b";
    c.shadowBlur = 5;
    c.font = "bold 9px 'JetBrains Mono', monospace";
    c.textAlign = "center";
    c.fillText("PPi", ppi.x, ppi.y);
    c.restore();
  }

  c.restore();
}

// Molecular Ambient Particles - Dynamically redistributed on layout/fullscreen changes
function updateAndDrawParticles(w, h, dt) {
  const c = ctx;
  if (!state.particles) state.particles = [];

  // Track canvas width and height changes (e.g. sidebar toggle or fullscreen)
  if (!state._lastParticleW) state._lastParticleW = w;
  if (!state._lastParticleH) state._lastParticleH = h;

  if (Math.abs(w - state._lastParticleW) > 10) {
    const scaleX = w / Math.max(state._lastParticleW, 1);
    for (let p of state.particles) {
      p.x *= scaleX;
      if (p.x < 0 || p.x > w) p.x = Math.random() * w;
    }
    state._lastParticleW = w;
  }
  if (Math.abs(h - state._lastParticleH) > 10) {
    const scaleY = h / Math.max(state._lastParticleH, 1);
    for (let p of state.particles) {
      p.y *= scaleY;
      if (p.y < 0 || p.y > h) p.y = Math.random() * h;
    }
    state._lastParticleH = h;
  }

  // Density-adjusted particle count based on canvas width
  const targetCount = Math.max(28, Math.floor(w / 38));
  while (state.particles.length < targetCount) {
    state.particles.push({
      x: Math.random() * w,
      y: Math.random() * h,
      vx: (Math.random() - 0.5) * 22,
      vy: (Math.random() - 0.5) * 22,
      size: Math.random() * 2 + 1,
      color: Math.random() > 0.5 ? "rgba(6, 182, 212, 0.45)" : "rgba(168, 85, 247, 0.45)"
    });
  }
  if (state.particles.length > targetCount + 10) {
    state.particles.length = targetCount;
  }

  c.save();
  for (let p of state.particles) {
    p.x += p.vx * dt;
    p.y += p.vy * dt;
    if (p.x < 0) p.x = w;
    if (p.x > w) p.x = 0;
    if (p.y < 0) p.y = h;
    if (p.y > h) p.y = 0;

    c.fillStyle = p.color;
    c.beginPath();
    c.arc(p.x, p.y, p.size, 0, Math.PI * 2);
    c.fill();
  }
  c.restore();
}
