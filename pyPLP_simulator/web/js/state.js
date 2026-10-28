/**
 * pyPLP Simulator - Simulation State, Viewport State & Settings
 */

// --- 2. MOLECULAR SIMULATION STATE ---
let state = {
  target: "",
  plp: "",
  primer: "",
  probe: "",
  temp: 37.0,
  targetType: "DNA",
  
  leftArmLen: 15,
  rightArmLen: 15,
  leftArmSeq: "",
  rightArmSeq: "",
  backboneSeq: "",
  
  leftTargetIdx: -1,
  rightTargetIdx: -1,
  nickGap: 0,
  isLigatable: false,
  hasNickMismatch: false,
  leftTm: 0,
  rightTm: 0,
  deltaTm: 0,
  
  // RCA Primer Binding State
  isPrimerComplementary: false,
  primerMatchPct: 0,
  primerDockIdx: -1,
  primerTm: 0,
  
  // Detection Probe Binding State
  isProbeComplementary: false,
  probeMatchPct: 0,
  probeDockIdx: -1,
  probeTm: 0,
  probeRepeatSites: 0,
  
  currentStage: 1,
  stageProgress: 0.0,
  isPlaying: false,
  playbackSpeed: 1.0,
  
  rcaRotationAngle: 0,
  concatemerRepeats: 0,
  probesBound: 0,
  fluoroSignal: 0,
  fluoroHistory: [0],
  
  // Interactive Protocol Pipeline State (Initially Empty Buffer)
  protocol: {
    hasTarget: false,
    hasPlp: false,
    hasLigase: false,
    hasPrimer: false,
    isPrimerAnnealed: false,
    hasPhi29: false,
    hasDntp: false,
    targetDisplaced: false,
    hasProbe: false,
    isProbeBound: false,
    userRepeats: 3
  },

  dntpPool: [],
  releasedPpi: [],
  particles: [],
  showDntps: false,
  stage5CircleRadius: 260
};

function toggleShowDntps(checked) {
  state.showDntps = !!checked;
  const cb = document.getElementById('toggleDntpCheckbox');
  if (cb) cb.checked = state.showDntps;
}
window.toggleShowDntps = toggleShowDntps;

function setCircleRadius(val) {
  let r = parseInt(val);
  if (isNaN(r)) r = 130;
  r = Math.max(50, Math.min(260, r));
  state.stage5CircleRadius = r;
  
  const inputEl = document.getElementById('input-circle-radius');
  if (inputEl) inputEl.value = r;
  
  const sliderEl = document.getElementById('slider-circle-radius');
  if (sliderEl) sliderEl.value = r;
  
  const valEl = document.getElementById('circle-radius-val');
  if (valEl) valEl.innerText = `${r} px`;
}
window.setCircleRadius = setCircleRadius;

function changeCircleRadius(delta) {
  const current = (state.stage5CircleRadius != null) ? state.stage5CircleRadius : 260;
  setCircleRadius(current + delta);
}
window.changeCircleRadius = changeCircleRadius;

// --- 3. PAN & ZOOM VIEWPORT STATE ---
let viewState = {
  zoom: 1.0,
  panX: 0,
  panY: 0,
  isDragging: false,
  dragStartX: 0,
  dragStartY: 0
};

function zoomCanvas(factor, centerPoint = null) {
  const targetCanvas = (typeof canvas !== 'undefined' && canvas) ? canvas : document.getElementById('reactionCanvas');
  const prevZoom = viewState.zoom;
  let newZoom = Math.min(Math.max(viewState.zoom * factor, 0.5), 6.0);
  
  const width = targetCanvas ? (targetCanvas.width / (window.devicePixelRatio || 1)) : 800;
  const height = targetCanvas ? (targetCanvas.height / (window.devicePixelRatio || 1)) : 600;
  const cx = centerPoint ? centerPoint.x : width * 0.5;
  const cy = centerPoint ? centerPoint.y : height * 0.5;

  viewState.panX = cx - (cx - viewState.panX) * (newZoom / prevZoom);
  viewState.panY = cy - (cy - viewState.panY) * (newZoom / prevZoom);
  viewState.zoom = newZoom;
  updateZoomUI();
}

function resetZoomPan() {
  viewState.zoom = 1.0;
  viewState.panX = 0;
  viewState.panY = 0;
  updateZoomUI();
}

function updateZoomUI() {
  const badge = document.getElementById('zoom-badge');
  if (badge) {
    const pct = Math.round(viewState.zoom * 100);
    badge.innerText = `${pct}%`;
  }
  const hint = document.getElementById('zoom-hint');
  if (hint) {
    const pct = Math.round(viewState.zoom * 100);
    if (shouldDrawBases()) {
      hint.innerHTML = `<span class="text-emerald-400 font-bold">● Bases Visible (${pct}%)</span> • Spacing: ${renderSettings.spacing}px`;
    } else {
      hint.innerHTML = `Scroll: Zoom • Drag: Pan • F: Fullscreen`;
    }
  }
}

// --- VISUALIZATION & NUCLEOTIDE SPACING SETTINGS ---
let renderSettings = {
  spacing: 32, // Distance between consecutive nucleotides in px
  radius: 6,  // Nucleotide circle radius in px
  lettersMode: 'auto' // 'auto' (zoom >= 1.25) or 'always'
};

function shouldDrawBases() {
  if (renderSettings.lettersMode === 'always') return true;
  return viewState.zoom >= 1.25;
}

function toggleSettingsModal() {
  const modal = document.getElementById('settings-modal');
  if (modal) {
    modal.classList.toggle('hidden');
  }
}

function updateSettings(key, value) {
  if (key === 'spacing') {
    renderSettings.spacing = parseInt(value);
    const spVal = document.getElementById('spacing-val');
    if (spVal) spVal.innerText = `${value} px`;
  } else if (key === 'radius') {
    renderSettings.radius = parseInt(value);
    const radVal = document.getElementById('radius-val');
    if (radVal) radVal.innerText = `${value} px`;
  } else if (key === 'lettersMode') {
    renderSettings.lettersMode = value;
    const btnAuto = document.getElementById('btn-letters-auto');
    const btnAlways = document.getElementById('btn-letters-always');
    if (value === 'auto') {
      if (btnAuto) btnAuto.className = "px-2 py-1 rounded bg-cyan-600 text-white font-semibold text-[11px] transition";
      if (btnAlways) btnAlways.className = "px-2 py-1 rounded bg-slate-800 text-slate-400 hover:text-white font-semibold text-[11px] transition";
    } else {
      if (btnAuto) btnAuto.className = "px-2 py-1 rounded bg-slate-800 text-slate-400 hover:text-white font-semibold text-[11px] transition";
      if (btnAlways) btnAlways.className = "px-2 py-1 rounded bg-cyan-600 text-white font-semibold text-[11px] transition";
    }
  }
  updateZoomUI();
}

function applySpacingPreset(spacing, radius) {
  const sSlider = document.getElementById('slider-spacing');
  const rSlider = document.getElementById('slider-radius');
  if (sSlider) sSlider.value = spacing;
  if (rSlider) rSlider.value = radius;
  updateSettings('spacing', spacing);
  updateSettings('radius', radius);
}
