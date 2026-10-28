/**
 * pyPLP Simulator - Stage Management & HUD Controls
 */

const STAGE_INFO = {
  1: {
    num: "STAGE 1",
    title: "Initial Target & Padlock Probe Approaching",
    desc: "Target strand (DNA/RNA) and linear Padlock Probe with 5' Left Arm and 3' Right Arm freely floating in solution."
  },
  2: {
    num: "STAGE 2",
    title: "Target Hybridization & Horseshoe Loop Formation",
    desc: "Left and Right arms bind juxtaposed to the target strand, forming an open horseshoe loop bringing the 3'-OH and 5'-PO4 into contiguous alignment."
  },
  3: {
    num: "STAGE 3",
    title: "Enzymatic Nick Interrogation & Circularization",
    desc: "Ligase (T4 DNA Ligase or SplintR) interrogates the junction. If a 0-bp flush nick is present, covalent phosphodiester bond seals the padlock into a closed circle."
  },
  4: {
    num: "STAGE 4",
    title: "RCA Amplification Primer Annealing",
    desc: "Sequence-specific RCA primer anneals to its complementary docking cassette within the circularized padlock probe backbone."
  },
  5: {
    num: "STAGE 5",
    title: "Phi29 Polymerase Rolling Circle Amplification (RCA)",
    desc: "Phi29 DNA Polymerase engages the primer 3'-end and continuously rolls around the closed circle, spooling out a long tandem concatemer ssDNA."
  },
  6: {
    num: "STAGE 6",
    title: "Detection Probe Hybridization & Fluorescence Emission",
    desc: "5'-FAM fluorescent detection probes hybridize periodically to each repeated cassette along the concatemer, emitting intense fluorescent reporter signal."
  }
};

function updateHUD() {
  const targetLenEl = document.getElementById('target-len');
  if (targetLenEl) targetLenEl.innerText = `${state.target.length} nt`;
  const plpLenEl = document.getElementById('plp-len');
  if (plpLenEl) plpLenEl.innerText = `${state.plp.length} nt`;
  const leftArmTagEl = document.getElementById('left-arm-tag');
  if (leftArmTagEl) leftArmTagEl.innerText = `${state.leftArmSeq.length} nt`;
  const rightArmTagEl = document.getElementById('right-arm-tag');
  if (rightArmTagEl) rightArmTagEl.innerText = `${state.rightArmSeq.length} nt`;
  const backboneTagEl = document.getElementById('backbone-tag');
  if (backboneTagEl) backboneTagEl.innerText = `${state.backboneSeq.length} nt`;
  const primerLenEl = document.getElementById('primer-len');
  if (primerLenEl) primerLenEl.innerText = `${state.primer.length} nt`;
  const probeLenEl = document.getElementById('probe-len');
  if (probeLenEl) probeLenEl.innerText = `${state.probe.length} nt`;

  const nickElem = document.getElementById('stat-nick');
  const hudNick = document.getElementById('hud-nick-status');
  if (nickElem && hudNick) {
    if (state.isLigatable) {
      nickElem.className = "text-sm font-bold text-emerald-400 font-mono mt-0.5";
      nickElem.innerText = "0 bp Flush Nick";
      hudNick.className = "px-2 py-0.5 rounded bg-emerald-500/20 text-emerald-400 font-mono text-xs font-bold border border-emerald-500/30";
      hudNick.innerText = "0-BP FLUSH NICK";
    } else if (state.hasNickMismatch) {
      nickElem.className = "text-sm font-bold text-rose-400 font-mono mt-0.5";
      nickElem.innerText = "3' Nick Mismatch (Arrested)";
      hudNick.className = "px-2 py-0.5 rounded bg-rose-500/20 text-rose-400 font-mono text-xs font-bold border border-rose-500/30";
      hudNick.innerText = "3' MISMATCH";
    } else {
      nickElem.className = "text-sm font-bold text-rose-400 font-mono mt-0.5";
      nickElem.innerText = state.nickGap === 999 ? "No Alignment" : `${state.nickGap} bp Gap (Non-ligatable)`;
      hudNick.className = "px-2 py-0.5 rounded bg-rose-500/20 text-rose-400 font-mono text-xs font-bold border border-rose-500/30";
      hudNick.innerText = state.nickGap === 999 ? "NO ALIGNMENT" : `${state.nickGap}-BP GAP`;
    }
  }

  const statTm = document.getElementById('stat-tm');
  if (statTm) statTm.innerText = `${state.leftTm}°C / ${state.rightTm}°C`;
  const statDeltaTm = document.getElementById('stat-deltatm');
  if (statDeltaTm) statDeltaTm.innerText = `ΔTm: ${state.deltaTm}°C (≤ 6°C)`;

  // RCA Primer Diagnostic Card
  const primerElem = document.getElementById('stat-primer');
  const primerTmElem = document.getElementById('stat-primer-tm');
  if (primerElem && primerTmElem) {
    if (state.isPrimerComplementary) {
      primerElem.className = "text-sm font-bold text-emerald-400 font-mono mt-0.5";
      primerElem.innerText = `Bound (${state.primerMatchPct}% match)`;
      primerTmElem.innerText = `Tm: ${state.primerTm}°C (Annealing Viable)`;
    } else {
      primerElem.className = "text-sm font-bold text-rose-400 font-mono mt-0.5";
      primerElem.innerText = `Mismatch (${state.primerMatchPct}% match)`;
      primerTmElem.innerText = `No Binding Site (RCA Blocked)`;
    }
  }
  
  // Overall Viability Diagnostic Card
  const viabilityElem = document.getElementById('stat-viability');
  const viabilityDesc = document.getElementById('stat-viability-desc');
  const armsViable = (state.leftTm >= state.temp && state.rightTm >= state.temp);

  if (viabilityElem && viabilityDesc) {
    if (state.isLigatable && armsViable && state.isPrimerComplementary && state.isProbeComplementary) {
      viabilityElem.className = "text-sm font-bold text-emerald-400 font-mono mt-0.5 flex items-center gap-1";
      viabilityElem.innerHTML = `<svg class="w-4 h-4 text-emerald-400" fill="currentColor" viewBox="0 0 20 20"><path fill-rule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clip-rule="evenodd"/></svg> PASS (Full Reaction Viable)`;
      viabilityDesc.innerText = state.targetType === 'RNA' ? "SplintR + Phi29 + Probe Validated" : "T4 + Phi29 + Probe Validated";
    } else {
      viabilityElem.className = "text-sm font-bold text-rose-400 font-mono mt-0.5 flex items-center gap-1";
      viabilityElem.innerHTML = `<svg class="w-4 h-4 text-rose-400" fill="currentColor" viewBox="0 0 20 20"><path fill-rule="evenodd" d="M10 18a8 8 0 100-16 8 8 0 000 16zm3.707-9.293a1 1 0 00-1.414-1.414L9 10.586 7.707 9.293a1 1 0 00-1.414 1.414l2 2a1 1 0 001.414 0l4-4z" clip-rule="evenodd"/></svg> FAIL (Reaction Arrested)`;
      if (!state.isLigatable) viabilityDesc.innerText = "Stage 3 Ligation Gap / Mismatch";
      else if (!state.isPrimerComplementary) viabilityDesc.innerText = "Stage 4 RCA Primer Non-Complementary";
      else if (!state.isProbeComplementary) viabilityDesc.innerText = "Stage 6 Probe Non-Complementary (0 a.u.)";
      else viabilityDesc.innerText = "Arm Tm below reaction temperature";
    }
  }

  const hudReactionTemp = document.getElementById('hud-reaction-temp');
  if (hudReactionTemp) hudReactionTemp.innerText = `${state.temp}°C`;
}

function setStage(s) {
  state.currentStage = s;
  state.stageProgress = 0.0;
  if (s >= 2) {
    state.protocol.hasTarget = true;
    state.protocol.hasPlp = true;
  }
  if (s >= 3) {
    state.protocol.hasLigase = true;
  }
  updateProtocolUI();
  updateStageButtons();
  updateStageText();
  
  const fluoroCard = document.getElementById('fluoro-card');
  if (fluoroCard) {
    if (s >= 5) {
      fluoroCard.classList.remove('hidden');
    } else {
      fluoroCard.classList.add('hidden');
    }
  }

  const radContainer = document.getElementById('ctrl-circle-radius-container');
  if (radContainer) {
    if (s >= 5) {
      radContainer.classList.add('border-cyan-500/80', 'shadow-[0_0_10px_rgba(6,182,212,0.3)]', 'bg-cyan-950/30');
      radContainer.classList.remove('border-slate-700', 'bg-slate-800/90');
    } else {
      radContainer.classList.remove('border-cyan-500/80', 'shadow-[0_0_10px_rgba(6,182,212,0.3)]', 'bg-cyan-950/30');
      radContainer.classList.add('border-slate-700', 'bg-slate-800/90');
    }
  }
}

function nextStage() {
  if (state.currentStage < 6) {
    setStage(state.currentStage + 1);
  }
}

function prevStage() {
  if (state.currentStage > 1) {
    setStage(state.currentStage - 1);
  }
}

function updateStageButtons() {
  for (let i = 1; i <= 6; i++) {
    const btn = document.getElementById(`btn-stage-${i}`);
    if (!btn) continue;
    if (i === state.currentStage) {
      btn.className = "stage-btn px-3 py-1.5 rounded-lg text-xs font-semibold badge-stage-active transition flex items-center gap-1.5";
    } else if (i < state.currentStage) {
      btn.className = "stage-btn px-3 py-1.5 rounded-lg text-xs font-semibold text-emerald-400 bg-emerald-950/40 border border-emerald-800/40 hover:bg-slate-800 transition flex items-center gap-1.5";
    } else {
      btn.className = "stage-btn px-3 py-1.5 rounded-lg text-xs font-semibold text-slate-400 hover:bg-slate-800 transition flex items-center gap-1.5";
    }
  }
}

function updateStageText() {
  const stageConfig = (activeTextConfig && activeTextConfig.stages && activeTextConfig.stages[state.currentStage]) || {};
  const fallbackInfo = STAGE_INFO[state.currentStage] || {};
  let num = fallbackInfo.num || `STAGE ${state.currentStage}`;
  let title = stageConfig.title !== undefined ? stageConfig.title : (fallbackInfo.title || "");
  let desc = stageConfig.desc !== undefined ? stageConfig.desc : (fallbackInfo.desc || "");

  if (state.currentStage === 3 && !state.isLigatable) {
    num = "STAGE 3 (ARRESTED)";
    title = state.hasNickMismatch 
      ? getCanvasText('ligation', 'arrestedMismatchBanner', "LIGATION ARRESTED: 3' MISMATCH DETECTED") 
      : (state.nickGap !== 999 ? getCanvasText('ligation', 'arrestedGapBanner', "LIGATION ARRESTED: {gap}-BP GAP DETECTED", { gap: state.nickGap }) : getCanvasText('ligation', 'arrestedNoAlignBanner', "LIGATION ARRESTED: ARMS DO NOT ALIGN"));
    desc = "Ligase cannot bridge the non-flush junction or mismatch. Padlock probe remains linear; circularization is blocked.";
  } else if (state.currentStage === 4 && !state.isPrimerComplementary) {
    num = "STAGE 4 (ARRESTED)";
    title = getCanvasText('rcaPrimer', 'primerArrestedHeader', "PRIMER ANNEALING ARRESTED: NO COMPLEMENTARY BINDING SITE");
    desc = getCanvasText('rcaPrimer', 'primerArrestedDetail', "Primer sequence matches 0 cassettes on circular template ({matchPct}% match < 70% required)", { matchPct: state.primerMatchPct });
  } else if (state.currentStage === 5 && (!state.isLigatable || !state.isPrimerComplementary)) {
    num = "STAGE 5 (BLOCKED)";
    const reason = !state.isLigatable ? "Ligation Failed" : "Primer Mismatch";
    title = getCanvasText('rollingCircle', 'rcaBlockedHeader', "RCA INITIATION BLOCKED: {reason}", { reason: reason });
    desc = !state.isLigatable 
      ? getCanvasText('rollingCircle', 'rcaReasonLigation', "Rolling circle amplification requires a closed circular template. Ligation was arrested.")
      : getCanvasText('rollingCircle', 'rcaReasonPrimer', "Phi29 Polymerase requires an annealed primer 3'-OH terminus to initiate.");
  } else if (state.currentStage === 6 && !state.isProbeComplementary) {
    num = "STAGE 6 (NO SIGNAL)";
    title = getCanvasText('detection', 'probeFailedHeader', "PROBE HYBRIDIZATION FAILED: PROBE NOT COMPLEMENTARY TO CONCATEMER");
    desc = getCanvasText('detection', 'probeFailedDetail', "Probe sequence has 0 complementary target cassettes ({matchPct}% match < 70%); 0 a.u. signal", { matchPct: state.probeMatchPct });
  }

  if (activeTextConfig && activeTextConfig._hidden && activeTextConfig._hidden[`stages.${state.currentStage}`]) {
    title = "";
    desc = "";
  }

  const hudStageNum = document.getElementById('hud-stage-num');
  if (hudStageNum) hudStageNum.innerText = num;
  const hudTitle = document.getElementById('hud-title');
  if (hudTitle) hudTitle.innerText = title;
  const hudDesc = document.getElementById('hud-desc');
  if (hudDesc) hudDesc.innerText = desc;
}

function updatePlayPauseButton() {
  const playIcon = document.getElementById('icon-play');
  const pauseIcon = document.getElementById('icon-pause');
  const playText = document.getElementById('btn-play-text');
  if (!playIcon || !pauseIcon || !playText) return;
  if (state.isPlaying) {
    playIcon.classList.add('hidden');
    pauseIcon.classList.remove('hidden');
    playText.innerText = "Pause";
  } else {
    playIcon.classList.remove('hidden');
    pauseIcon.classList.add('hidden');
    playText.innerText = "Play";
  }
}

function togglePlay() {
  state.isPlaying = !state.isPlaying;
  updatePlayPauseButton();
}

function setSpeed(v) {
  state.playbackSpeed = parseFloat(v);
  const slider = document.getElementById('speedSlider');
  if (slider) slider.value = v;
  const valSpan = document.getElementById('speedVal');
  if (valSpan) valSpan.innerText = `${parseFloat(v).toFixed(2)}x`;
}

function resetSimulation() {
  state.currentStage = 1;
  state.stageProgress = 0.0;
  state.rcaRotationAngle = 0;
  state.concatemerRepeats = 0;
  state.fluoroSignal = 0;
  state.fluoroHistory = [0];
  state.particles = [];
  updateStageButtons();
  updateStageText();
  const fluoroCard = document.getElementById('fluoro-card');
  if (fluoroCard) fluoroCard.classList.add('hidden');
  updateProtocolUI();
}
