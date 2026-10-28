/**
 * pyPLP Simulator - Application Entrypoint & Event Handling
 */

// Canvas & Rendering Context Setup
const canvas = document.getElementById('reactionCanvas');
const ctx = canvas.getContext('2d');
const fCanvas = document.getElementById('fluoroCanvas');
const fCtx = fCanvas ? fCanvas.getContext('2d') : null;

function resizeCanvas() {
  if (!canvas || !canvas.parentElement) return;
  const rect = canvas.parentElement.getBoundingClientRect();
  const dpr = window.devicePixelRatio || 1;
  const newW = Math.round(rect.width * dpr);
  const newH = Math.round(rect.height * dpr);
  if (canvas.width !== newW || canvas.height !== newH) {
    canvas.width = newW;
    canvas.height = newH;
  }
}

window.addEventListener('resize', resizeCanvas);

// Automatically track layout shifts, sidebar collapse/expand, and responsive resizing
if (window.ResizeObserver && canvas && canvas.parentElement) {
  const resizeObserver = new ResizeObserver(() => {
    resizeCanvas();
  });
  resizeObserver.observe(canvas.parentElement);
}

// Mouse Wheel Zoom centered on cursor
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  const rect = canvas.getBoundingClientRect();
  const mouseX = e.clientX - rect.left;
  const mouseY = e.clientY - rect.top;
  const factor = e.deltaY < 0 ? 1.15 : 0.87;
  zoomCanvas(factor, { x: mouseX, y: mouseY });
}, { passive: false });

// Drag to Pan
canvas.addEventListener('mousedown', (e) => {
  viewState.isDragging = true;
  viewState.dragStartX = e.clientX - viewState.panX;
  viewState.dragStartY = e.clientY - viewState.panY;
  canvas.classList.add('cursor-grabbing');
  canvas.classList.remove('cursor-grab');
});

window.addEventListener('mousemove', (e) => {
  if (viewState.isDragging) {
    viewState.panX = e.clientX - viewState.dragStartX;
    viewState.panY = e.clientY - viewState.dragStartY;
  }
});

window.addEventListener('mouseup', () => {
  if (viewState.isDragging) {
    viewState.isDragging = false;
    canvas.classList.remove('cursor-grabbing');
    canvas.classList.add('cursor-grab');
  }
});

// Touch Support for Mobile / Tablets
let initialTouchDist = null;
canvas.addEventListener('touchstart', (e) => {
  if (e.touches.length === 1) {
    viewState.isDragging = true;
    viewState.dragStartX = e.touches[0].clientX - viewState.panX;
    viewState.dragStartY = e.touches[0].clientY - viewState.panY;
  } else if (e.touches.length === 2) {
    initialTouchDist = Math.hypot(
      e.touches[0].clientX - e.touches[1].clientX,
      e.touches[0].clientY - e.touches[1].clientY
    );
  }
});

canvas.addEventListener('touchmove', (e) => {
  if (e.touches.length === 1 && viewState.isDragging) {
    viewState.panX = e.touches[0].clientX - viewState.dragStartX;
    viewState.panY = e.touches[0].clientY - viewState.dragStartY;
  } else if (e.touches.length === 2 && initialTouchDist) {
    const currentDist = Math.hypot(
      e.touches[0].clientX - e.touches[1].clientX,
      e.touches[0].clientY - e.touches[1].clientY
    );
    const factor = currentDist / initialTouchDist;
    zoomCanvas(factor > 1 ? 1.05 : 0.95);
    initialTouchDist = currentDist;
  }
});

canvas.addEventListener('touchend', () => {
  viewState.isDragging = false;
  initialTouchDist = null;
});

// --- FULLSCREEN CONTROLS & SHORTCUTS ---
function toggleFullscreen() {
  const elem = document.getElementById('canvas-container');
  if (!elem) return;
  if (!document.fullscreenElement && !document.webkitFullscreenElement) {
    if (elem.requestFullscreen) {
      elem.requestFullscreen();
    } else if (elem.webkitRequestFullscreen) {
      elem.webkitRequestFullscreen();
    }
  } else {
    if (document.exitFullscreen) {
      document.exitFullscreen();
    } else if (document.webkitExitFullscreen) {
      document.webkitExitFullscreen();
    }
  }
}

function handleFullscreenChange() {
  const isFS = !!(document.fullscreenElement || document.webkitFullscreenElement);
  const maxIcon = document.getElementById('icon-maximize');
  const minIcon = document.getElementById('icon-minimize');
  if (maxIcon && minIcon) {
    if (isFS) {
      maxIcon.classList.add('hidden');
      minIcon.classList.remove('hidden');
    } else {
      maxIcon.classList.remove('hidden');
      minIcon.classList.add('hidden');
    }
  }
  setTimeout(resizeCanvas, 80);
}

document.addEventListener('fullscreenchange', handleFullscreenChange);
document.addEventListener('webkitfullscreenchange', handleFullscreenChange);

// Global Keyboard Shortcuts
window.addEventListener('keydown', (e) => {
  const tag = (document.activeElement && document.activeElement.tagName) || '';
  if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;

  if (e.key === 'f' || e.key === 'F') {
    e.preventDefault();
    toggleFullscreen();
  } else if (e.key === ' ') {
    e.preventDefault();
    togglePlay();
  } else if (e.key === 'ArrowRight') {
    e.preventDefault();
    nextStage();
  } else if (e.key === 'ArrowLeft') {
    e.preventDefault();
    prevStage();
  } else if (e.key === '0') {
    e.preventDefault();
    resetZoomPan();
  } else if (e.key === '+' || e.key === '=') {
    e.preventDefault();
    zoomCanvas(1.25);
  } else if (e.key === '-' || e.key === '_') {
    e.preventDefault();
    zoomCanvas(0.8);
  }
});

// --- ANIMATION LOOP & STAGE DRAWERS ---
let lastTime = performance.now();

function renderLoop(now) {
  const dt = (now - lastTime) / 1000;
  lastTime = now;

  const dpr = window.devicePixelRatio || 1;
  const width = canvas.width / dpr;
  const height = canvas.height / dpr;

  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);

  // Apply Pan & Zoom Viewport Transformation
  ctx.save();
  ctx.translate(viewState.panX, viewState.panY);
  ctx.scale(viewState.zoom, viewState.zoom);

  drawBackgroundGrid(width, height);

  if (state.isPlaying) {
    let rate = 0.45;
    if (state.currentStage === 5) {
      const uRepeats = Math.max(1, (state.protocol && state.protocol.userRepeats) || 3);
      if (state.stageProgress < 0.10) {
        // Phase 1: Phi29 docking onto primer 3'-OH (~2.0s)
        rate = 0.05;
      } else if (state.stageProgress < 0.42) {
        // Phase 2: Displacing target strand base-by-base (~3.5 bases/sec, ~14.5s for 51 bases)
        rate = 0.022;
      } else if (state.stageProgress < 0.56) {
        // Phase 3: Target displaced & template morphing to circle (~2.8s)
        rate = 0.05;
      } else {
        // Phase 4: Continuous rolling circle amplification (matched to ~3.5 bases/sec)
        rate = 0.021 / uRepeats;
      }
    }
    state.stageProgress += dt * rate * state.playbackSpeed;
    if (state.stageProgress >= 1.0) {
      if (state.currentStage === 3) {
        // Pause when ligation completes; do NOT auto-add primer or auto-advance
        state.stageProgress = 1.0;
        state.isPlaying = false;
        updatePlayPauseButton();
      } else if (state.currentStage === 4) {
        // Pause at end of Stage 4 — user clicks + Add phi29 to proceed
        state.stageProgress = 1.0;
        state.isPlaying = false;
        state.protocol.isPrimerAnnealed = true;
        updatePlayPauseButton();
        updateProtocolUI();
      } else if (state.currentStage === 5) {
        // RCA Repeats completed — pause playback and stop rotation
        state.stageProgress = 1.0;
        state.isPlaying = false;
        updatePlayPauseButton();
      } else if (state.currentStage < 6) {
        setStage(state.currentStage + 1);
      } else {
        state.stageProgress = 1.0;
        state.isPlaying = false;
        updatePlayPauseButton();
      }
    }
  }

  switch (state.currentStage) {
    case 1:
      drawStage1(width, height, state.stageProgress);
      break;
    case 2:
      drawStage2(width, height, state.stageProgress);
      break;
    case 3:
      drawStage3(width, height, state.stageProgress);
      break;
    case 4:
      drawStage4(width, height, state.stageProgress);
      break;
    case 5:
      drawStage5(width, height, state.stageProgress, dt);
      break;
    case 6:
      drawStage6(width, height, state.stageProgress, dt);
      break;
  }

  ctx.restore(); // Restore pan/zoom

  // Ambient floating particles rendered in screen-space across the full canvas viewport
  updateAndDrawParticles(width, height, dt);

  if (state.currentStage >= 5 && fCanvas && fCtx) {
    drawFluorometerGraph();
  }

  requestAnimationFrame(renderLoop);
}

// --- USER INPUTS & PRESET LOADER ---
function loadPreset(key) {
  const p = PRESETS[key];
  if (!p) return;

  const inTarget = document.getElementById('input-target');
  const inPlp = document.getElementById('input-plp');
  const inPrimer = document.getElementById('input-primer');
  const inProbe = document.getElementById('input-probe');
  const inTemp = document.getElementById('input-temp');
  const inType = document.getElementById('select-target-type');

  if (inTarget) inTarget.value = p.target;
  if (inPlp) inPlp.value = p.plp;
  if (inPrimer) inPrimer.value = p.primer;
  if (inProbe) inProbe.value = p.probe;
  if (inTemp) inTemp.value = p.temp;
  if (inType) inType.value = p.type;

  state.protocol = {
    hasTarget: true,
    hasPlp: true,
    hasLigase: true,
    hasPrimer: true,
    isPrimerAnnealed: false,
    hasPhi29: true,
    targetDisplaced: false,
    userRepeats: (state.protocol && state.protocol.userRepeats) || 3
  };
  const repInput = document.getElementById('input-rca-repeats');
  if (repInput) repInput.value = state.protocol.userRepeats;

  applyUserInputs();
  resetSimulation();
  setStage(1);
}

function toggleReactionParams() {
  const body = document.getElementById('rxn-params-body');
  const chevron = document.getElementById('rxn-params-chevron');
  const statusBadge = document.getElementById('rxn-params-status');
  if (!body) return;
  const isHidden = body.classList.contains('hidden');
  if (isHidden) {
    body.classList.remove('hidden');
    if (chevron) chevron.style.transform = 'rotate(0deg)';
    if (statusBadge) statusBadge.innerText = 'Interactive Inputs';
  } else {
    body.classList.add('hidden');
    if (chevron) chevron.style.transform = 'rotate(-90deg)';
    if (statusBadge) statusBadge.innerText = 'Click to Expand';
  }
}
window.toggleReactionParams = toggleReactionParams;

function toggleSidebar() {
  const sidebar = document.getElementById('sidebar-col');
  const mainCol = document.getElementById('main-content-col');
  const btn = document.getElementById('btn-toggle-sidebar');
  const btnText = document.getElementById('btn-sidebar-text');
  if (!sidebar || !mainCol) return;

  const isCollapsed = sidebar.classList.contains('hidden');
  if (isCollapsed) {
    // Expand sidebar: restore 8 cols + 4 cols
    sidebar.classList.remove('hidden');
    mainCol.classList.remove('lg:col-span-12');
    mainCol.classList.add('lg:col-span-8');
    if (btn) {
      btn.className = "px-3 py-1 text-xs font-semibold rounded-md bg-gradient-to-r from-cyan-600/40 to-blue-600/40 hover:from-cyan-600 hover:to-blue-600 border border-cyan-500/50 text-cyan-200 hover:text-white transition flex items-center gap-1.5 shadow";
    }
    if (btnText) btnText.innerText = "Parameters Panel";
  } else {
    // Collapse sidebar: mainCol becomes 12 cols, sidebar hidden
    sidebar.classList.add('hidden');
    mainCol.classList.remove('lg:col-span-8');
    mainCol.classList.add('lg:col-span-12');
    if (btn) {
      btn.className = "px-3 py-1 text-xs font-semibold rounded-md bg-slate-800 hover:bg-slate-700 border border-slate-700 text-slate-300 hover:text-white transition flex items-center gap-1.5 shadow";
    }
    if (btnText) btnText.innerText = "Show Parameters";
  }

  // Trigger canvas resize so the drawing fills the expanded width smoothly
  resizeCanvas();
  setTimeout(resizeCanvas, 50);
  setTimeout(resizeCanvas, 150);
  setTimeout(resizeCanvas, 320);
}
window.toggleSidebar = toggleSidebar;

function applyUserInputs() {
  const inTarget = document.getElementById('input-target');
  const inPlp = document.getElementById('input-plp');
  const inPrimer = document.getElementById('input-primer');
  const inProbe = document.getElementById('input-probe');
  const inTemp = document.getElementById('input-temp');
  const inType = document.getElementById('select-target-type');

  if (inTarget) state.target = inTarget.value.trim();
  if (inPlp) state.plp = inPlp.value.trim();
  if (inPrimer) state.primer = inPrimer.value.trim();
  if (inProbe) state.probe = inProbe.value.trim();
  if (inTemp) state.temp = parseFloat(inTemp.value);
  if (inType) state.targetType = inType.value;

  analyzeMolecules();
}

window.addEventListener('DOMContentLoaded', () => {
  resizeCanvas();
  loadPreset('mtb_dna');
  protocolResetBuffer();
  requestAnimationFrame(renderLoop);
});
