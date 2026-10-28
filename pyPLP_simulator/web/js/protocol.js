/**
 * pyPLP Simulator - Interactive Laboratory Protocol Actions & Controller
 */

function protocolAddTarget() {
  state.protocol.hasTarget = true;
  setStage(1);
  state.stageProgress = 0.0;
  updateProtocolUI();
}

function protocolAddPLP() {
  state.protocol.hasPlp = true;
  setStage(1);
  state.stageProgress = 0.0;
  updateProtocolUI();
}

function protocolAddLigase() {
  state.protocol.hasLigase = true;
  if (state.protocol.hasTarget && state.protocol.hasPlp) {
    if (state.currentStage === 2 && state.stageProgress >= 1.0) {
      setStage(3);
    } else {
      setStage(2);
    }
    state.stageProgress = 0.0;
    state.isPlaying = true;
    updatePlayPauseButton();
  }
  updateProtocolUI();
}

function protocolStartLigation() {
  state.protocol.hasTarget = true;
  state.protocol.hasPlp = true;
  state.protocol.hasLigase = true;
  setStage(2);
  state.stageProgress = 0.0;
  state.isPlaying = true;
  updatePlayPauseButton();
  updateProtocolUI();
}

function protocolAddPrimer() {
  state.protocol.hasTarget = true;
  state.protocol.hasPlp = true;
  state.protocol.hasLigase = true;
  state.protocol.hasPrimer = true;
  state.protocol.isPrimerAnnealed = false;
  setStage(4);
  state.stageProgress = 0.0;
  state.isPlaying = true;
  updatePlayPauseButton();
  updateProtocolUI();
}

function protocolAddPhi29() {
  // Absorbs primer annealing + phi29 activation in one click.
  state.protocol.hasTarget = true;
  state.protocol.hasPlp = true;
  state.protocol.hasLigase = true;
  state.protocol.hasPrimer = true;
  state.protocol.isPrimerAnnealed = true;
  state.protocol.hasPhi29 = true;
  state.protocol.hasDntp = true;
  setStage(5);
  state.stageProgress = 0.0;
  state.isPlaying = true;
  updatePlayPauseButton();
  updateProtocolUI();
}

function protocolAddProbe() {
  state.protocol.hasProbe = true;
  state.protocol.isProbeBound = true;
  setStage(6);
  state.stageProgress = 0.0;
  state.isPlaying = true;
  updatePlayPauseButton();
  updateProtocolUI();
}

function protocolBindProbe() {
  state.protocol.hasProbe = true;
  state.protocol.isProbeBound = true;
  updateProtocolUI();
}

function protocolResetBuffer() {
  state.protocol.hasTarget = false;
  state.protocol.hasPlp = false;
  state.protocol.hasLigase = false;
  state.protocol.hasPrimer = false;
  state.protocol.isPrimerAnnealed = false;
  state.protocol.hasPhi29 = false;
  state.protocol.hasDntp = false;
  state.protocol.targetDisplaced = false;
  state.protocol.hasProbe = false;
  state.protocol.isProbeBound = false;
  state.dntpPool = [];
  state.releasedPpi = [];
  state.isPlaying = false;
  setStage(1);
  state.stageProgress = 0.0;
  state.rcaRotationAngle = 0;
  state.concatemerRepeats = 0;
  state.fluoroSignal = 0;
  state.fluoroHistory = [0];
  updatePlayPauseButton();
  updateProtocolUI();
}

function setRcaRepeats(val) {
  let n = parseInt(val) || 3;
  n = Math.max(1, Math.min(20, n));
  state.protocol.userRepeats = n;
  const el = document.getElementById('input-rca-repeats');
  if (el) el.value = n;
}

function changeRcaRepeats(delta) {
  const current = state.protocol.userRepeats || 3;
  setRcaRepeats(current + delta);
}

function updateProtocolUI() {
  const p = state.protocol;
  const bTarget = document.getElementById('btn-proto-target');
  const bPlp = document.getElementById('btn-proto-plp');
  const bLigase = document.getElementById('btn-proto-ligase');
  const bPrimer = document.getElementById('btn-proto-primer');
  const bPhi29 = document.getElementById('btn-proto-phi29');
  const bDntp = document.getElementById('btn-proto-dntp');
  const bProbe = document.getElementById('btn-proto-probe');
  const bProbeBind = document.getElementById('btn-proto-probe-bind');

  if (bTarget) {
    bTarget.style.opacity = p.hasTarget ? "1.0" : "0.5";
    bTarget.innerHTML = p.hasTarget 
      ? `<span class="text-emerald-400 font-bold">✓</span> Target Added` 
      : `<span class="text-amber-400 font-bold">+</span> Add Target`;
  }
  if (bPlp) {
    bPlp.style.opacity = p.hasPlp ? "1.0" : "0.5";
    bPlp.innerHTML = p.hasPlp 
      ? `<span class="text-emerald-400 font-bold">✓</span> PLP Added` 
      : `<span class="text-purple-400 font-bold">+</span> Add PLP`;
  }
  if (bLigase) {
    bLigase.style.opacity = p.hasLigase ? "1.0" : "0.5";
    bLigase.innerHTML = p.hasLigase 
      ? `<span class="text-emerald-400 font-bold">✓</span> Ligase Added` 
      : `<span class="text-emerald-400 font-bold">+</span> Add Ligase`;
  }
  if (bPrimer) {
    bPrimer.style.opacity = p.hasPrimer ? "1.0" : "0.5";
    bPrimer.innerHTML = p.hasPrimer 
      ? (p.isPrimerAnnealed ? `<span class="text-emerald-400 font-bold">✓</span> Primer Annealed` : `<span class="text-rose-400 font-bold">●</span> Primer Added`)
      : `<span class="text-rose-400 font-bold">+</span> Add RCA Primer`;
  }
  if (bPhi29) {
    bPhi29.style.opacity = p.hasPhi29 ? "1.0" : "0.5";
    bPhi29.innerHTML = p.hasPhi29 
      ? `<span class="text-cyan-400 font-bold">✓</span> phi29 Active` 
      : `<span class="text-cyan-400 font-bold">+</span> Add phi29`;
  }
  if (bDntp) {
    bDntp.style.opacity = p.hasDntp ? "1.0" : "0.5";
    bDntp.innerHTML = p.hasDntp 
      ? `<span class="text-emerald-400 font-bold">✓</span> dNTPs Added` 
      : `<span class="text-emerald-400 font-bold">+</span> Add dNTPs`;
  }
  if (bProbe) {
    bProbe.style.opacity = p.hasProbe ? "1.0" : "0.5";
    bProbe.innerHTML = p.hasProbe 
      ? `<span class="text-yellow-400 font-bold">✓</span> Probe Added` 
      : `<span class="text-yellow-400 font-bold">+</span> Add Probe`;
  }
  if (bProbeBind) {
    bProbeBind.style.opacity = p.isProbeBound ? "1.0" : "0.5";
    bProbeBind.innerHTML = p.isProbeBound
      ? `<span class="text-emerald-400 font-bold">✓</span> Probes Bound`
      : `<span>⚡</span> Bind Probe`;
  }
}
