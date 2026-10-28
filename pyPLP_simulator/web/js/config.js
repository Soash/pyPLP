/**
 * pyPLP Simulator - Canvas Text Configuration Manager
 */

const DEFAULT_CANVAS_TEXT_CONFIG = {
  strands: {
    targetStrandLabel: "Target ({targetType}) 5' → 3'",
    leftArmLinearLabel: "5'-Left Arm",
    rightArmLinearLabel: "Right Arm-3'",
    backboneLinearLabel: "Backbone",
    leftArmHybridLabel: "Left Arm (5'-PO4)",
    rightArmHybridLabel: "Right Arm (3'-OH)",
    backboneLoopLabel: "Backbone Loop",
    nucleotide5PrimeTag: "5'",
    nucleotide3PrimeTag: "3'"
  },
  ligation: {
    nickFlushBadge: "0-bp Nick (3'-OH | 5'-PO4)",
    nickMismatchBadge: "3'-End Mismatch (Arrested)",
    nickGapBadge: "{gap}-bp Gap (Arrested)",
    nickNoAlignBadge: "No Alignment",
    enzymeDnaLigase: "Ligase",
    enzymeRnaLigase: "SplintR",
    sealedSuccessBanner: "COVALENT PHOSPHODIESTER BOND SEALED",
    arrestedMismatchBanner: "LIGATION ARRESTED: 3' MISMATCH DETECTED",
    arrestedGapBanner: "LIGATION ARRESTED: {gap}-BP GAP DETECTED",
    arrestedNoAlignBanner: "LIGATION ARRESTED: ARMS DO NOT ALIGN"
  },
  rcaPrimer: {
    circularPlpTitle: "Circular Padlock Probe (Closed Single Strand)",
    primerAnnealedBanner: "RCA PRIMER ANNEALED (3'-OH EXTENSION READY)",
    primerStrandLabel: "RCA Primer (5' → 3')",
    primerBlockedBadge: "✕ NON-COMPLEMENTARY PRIMER",
    primerArrestedHeader: "PRIMER ANNEALING ARRESTED: NO COMPLEMENTARY BINDING SITE",
    primerArrestedDetail: "Primer sequence matches 0 cassettes on circular template ({matchPct}% match < 70% required)"
  },
  rollingCircle: {
    phi29EnzymeName: "phi29",
    phi29InactiveSymbol: "✕",
    repeatBadgePrefix: "Repeat {repeatNum}",
    concatemerStrandLabel: "Single-Stranded Concatemer DNA (Antisense)",
    rcaBlockedHeader: "RCA INITIATION BLOCKED: {reason}",
    rcaReasonLigation: "Rolling circle amplification requires a closed circular template. Ligation was arrested.",
    rcaReasonPrimer: "Phi29 Polymerase requires an annealed primer 3'-OH terminus to initiate."
  },
  detection: {
    fluorophoreActiveBadge: "FAM",
    fluorophoreQuenchedBadge: "FAM (Off)",
    probeStrandLabel: "Detection Probe (5' → 3')",
    probeFailedHeader: "PROBE HYBRIDIZATION FAILED: PROBE NOT COMPLEMENTARY TO CONCATEMER",
    probeFailedDetail: "Probe sequence has 0 complementary target cassettes ({matchPct}% match < 70%); 0 a.u. signal",
    fluorometerGraphTitle: "REAL-TIME FLUOROMETER",
    fluorometerIntensityUnit: "a.u.",
    fluorometerTimeLabel: "t (s)"
  },
  stages: {
    1: {
      title: "Stage 1: Solution State (Target & Padlock Probe)",
      desc: "Single-stranded target sequence (DNA/RNA) and linear Padlock Probe with 5'-phosphate and 3'-hydroxyl termini float independently in hybridization buffer."
    },
    2: {
      title: "Stage 2: Target Recognition & Hybridization",
      desc: "Padlock probe arms recognize and anneal to contiguous target sites in an antiparallel orientation, looping the backbone and aligning termini into a nick junction."
    },
    3: {
      title: "Stage 3: Enzymatic Ligation & Topological Interlocking",
      desc: "DNA Ligase or SplintR inspects the nick. If perfect base-pairing is detected with zero gap, it catalyzes phosphodiester bond synthesis to yield a covalently closed circle."
    },
    4: {
      title: "Stage 4: RCA Primer Annealing",
      desc: "Sequence-specific RCA primer anneals to the complementary binding cassette on the circular padlock template, supplying the free 3'-OH terminus required for polymerase extension."
    },
    5: {
      title: "Stage 5: Continuous Phi29 Rolling Circle Amplification (RCA)",
      desc: "Phi29 DNA Polymerase engages the primed template and initiates isothermal rolling synthesis, displacing downstream strands to generate long tandem concatemers."
    },
    6: {
      title: "Stage 6: Concatemer Hybridization & Fluorescence Detection",
      desc: "Sequence-specific fluorescent detection probes anneal to repeating monomer cassettes along the concatemer, emitting measurable fluorescent signals."
    }
  }
};

// Deep merge helper
function deepMerge(target, source) {
  const output = Object.assign({}, target);
  if (source && typeof source === 'object') {
    Object.keys(source).forEach(key => {
      if (source[key] && typeof source[key] === 'object' && !Array.isArray(source[key])) {
        output[key] = deepMerge(target[key] || {}, source[key]);
      } else {
        output[key] = source[key];
      }
    });
  }
  return output;
}

let activeTextConfig = deepMerge(DEFAULT_CANVAS_TEXT_CONFIG, (typeof CANVAS_TEXT_CONFIG !== 'undefined') ? CANVAS_TEXT_CONFIG : {});

function getCanvasText(section, key, fallback, replacements = {}) {
  if (activeTextConfig && activeTextConfig._hidden && activeTextConfig._hidden[`${section}.${key}`]) {
    return "";
  }
  let val = fallback;
  if (activeTextConfig && activeTextConfig[section] && activeTextConfig[section][key] !== undefined) {
    val = activeTextConfig[section][key];
  }
  for (const [k, v] of Object.entries(replacements)) {
    val = String(val).replace(new RegExp(`\\{${k}\\}`, 'g'), v);
  }
  return val;
}

let currentCfgTab = 'strands';

function openTextConfigModal() {
  const modal = document.getElementById('text-config-modal');
  modal.classList.remove('hidden');
  renderConfigForm();
}

function closeTextConfigModal() {
  document.getElementById('text-config-modal').classList.add('hidden');
}

function switchConfigTab(tab) {
  currentCfgTab = tab;
  document.querySelectorAll('.cfg-tab').forEach(btn => {
    btn.className = "cfg-tab py-2.5 px-3 font-medium border-b-2 border-transparent text-slate-400 hover:text-slate-200 transition whitespace-nowrap";
  });
  const activeBtn = document.getElementById(`tab-btn-${tab}`);
  if (activeBtn) {
    activeBtn.className = "cfg-tab py-2.5 px-3 font-semibold border-b-2 border-indigo-500 text-indigo-400 transition whitespace-nowrap";
  }
  renderConfigForm();
}

const CONFIG_FIELD_META = {
  strands: [
    { key: "targetStrandLabel", label: "Target Strand Label", hint: "Supports {targetType}" },
    { key: "leftArmLinearLabel", label: "5' Left Arm Label (Stage 1 Linear)", hint: "Linear probe in solution" },
    { key: "rightArmLinearLabel", label: "3' Right Arm Label (Stage 1 Linear)", hint: "Linear probe in solution" },
    { key: "backboneLinearLabel", label: "Backbone Label (Stage 1 Linear)", hint: "Linear central backbone connector" },
    { key: "leftArmHybridLabel", label: "Left Arm Label (Stage 2/3 Hybridized)", hint: "Hybridized to left target arm" },
    { key: "rightArmHybridLabel", label: "Right Arm Label (Stage 2/3 Hybridized)", hint: "Hybridized to right target arm" },
    { key: "backboneLoopLabel", label: "Backbone Arched Loop Label", hint: "Arched over the hybridized duplex" },
    { key: "nucleotide5PrimeTag", label: "5' Terminus Badge", hint: "Drawn on first nucleotide node" },
    { key: "nucleotide3PrimeTag", label: "3' Terminus Badge", hint: "Drawn on last nucleotide node" }
  ],
  ligation: [
    { key: "nickFlushBadge", label: "0-bp Nick Badge", hint: "Indicator circle when flush nick aligns" },
    { key: "nickMismatchBadge", label: "3' Mismatch Badge", hint: "Indicator circle when -1 mismatch occurs" },
    { key: "nickGapBadge", label: "Nick Gap Badge", hint: "Supports {gap}" },
    { key: "nickNoAlignBadge", label: "No Alignment Badge", hint: "When target doesn't align" },
    { key: "enzymeDnaLigase", label: "DNA Ligase Enzyme Name", hint: "Drawn inside the enzyme circle" },
    { key: "enzymeRnaLigase", label: "RNA SplintR Ligase Enzyme Name", hint: "Drawn inside the enzyme circle" },
    { key: "sealedSuccessBanner", label: "Covalent Seal Success Banner", hint: "Drawn below ligated complex" },
    { key: "arrestedMismatchBanner", label: "Ligation Arrested (Mismatch Banner)", hint: "When 3' mismatch halts ligation" },
    { key: "arrestedGapBanner", label: "Ligation Arrested (Gap Banner)", hint: "Supports {gap}" },
    { key: "arrestedNoAlignBanner", label: "Ligation Arrested (No Alignment)", hint: "When arms cannot align" }
  ],
  rcaPrimer: [
    { key: "circularPlpTitle", label: "Circular Padlock Title", hint: "Header above closed circular template" },
    { key: "primerAnnealedBanner", label: "Primer Annealed Banner", hint: "Success banner when primer docks" },
    { key: "primerStrandLabel", label: "RCA Primer Strand Label", hint: "Label above annealed primer" },
    { key: "primerBlockedBadge", label: "Primer Blocked Warning Badge", hint: "Drawn over unannealed primer" },
    { key: "primerArrestedHeader", label: "Primer Arrested Header", hint: "Failure announcement title" },
    { key: "primerArrestedDetail", label: "Primer Arrested Detail", hint: "Supports {matchPct}" }
  ],
  rollingCircle: [
    { key: "phi29EnzymeName", label: "Phi29 Polymerase Name", hint: "Drawn inside active polymerase node" },
    { key: "phi29InactiveSymbol", label: "Inactive Polymerase Mark", hint: "Drawn inside arrested enzyme node" },
    { key: "repeatBadgePrefix", label: "Concatemer Repeat Marker", hint: "Supports {repeatNum}" },
    { key: "concatemerStrandLabel", label: "Concatemer Strand Label", hint: "Label above growing tail" },
    { key: "rcaBlockedHeader", label: "RCA Blocked Header", hint: "Supports {reason}" },
    { key: "rcaReasonLigation", label: "Arrest Reason (Ligation Failed)", hint: "Reason explanation" },
    { key: "rcaReasonPrimer", label: "Arrest Reason (Primer Failed)", hint: "Reason explanation" }
  ],
  detection: [
    { key: "fluorophoreActiveBadge", label: "Active Fluorophore Label", hint: "Drawn on pulsing fluorophore" },
    { key: "fluorophoreQuenchedBadge", label: "Quenched Fluorophore Label", hint: "Drawn on dark fluorophore" },
    { key: "probeStrandLabel", label: "Detection Probe Strand Label", hint: "Label for unbound probe" },
    { key: "probeFailedHeader", label: "Probe Hybridization Failed Title", hint: "When probe doesn't match" },
    { key: "probeFailedDetail", label: "Probe Hybridization Failed Detail", hint: "Supports {matchPct}" },
    { key: "fluorometerGraphTitle", label: "Fluorometer Graph Title", hint: "Mini graph card title" },
    { key: "fluorometerIntensityUnit", label: "Fluorometer Intensity Unit", hint: "Default: a.u." },
    { key: "fluorometerTimeLabel", label: "Fluorometer X-Axis Label", hint: "Default: t (s)" }
  ]
};

function renderConfigForm() {
  const container = document.getElementById('cfg-form-container');
  if (!container) return;

  if (currentCfgTab === 'stages') {
    let html = `
      <div class="flex items-center justify-between pb-2 mb-3 border-b border-slate-800">
        <span class="text-xs text-slate-400 font-medium">Stage Info Banners</span>
        <div class="flex items-center gap-2">
          <button type="button" onclick="setAllStageToggles(true)" class="px-2 py-0.5 text-[11px] rounded bg-slate-800 hover:bg-slate-700 text-cyan-300 transition">Show All</button>
          <button type="button" onclick="setAllStageToggles(false)" class="px-2 py-0.5 text-[11px] rounded bg-slate-800 hover:bg-slate-700 text-slate-400 hover:text-white transition">Hide All</button>
        </div>
      </div>
      <div class="space-y-4">
    `;
    for (let i = 1; i <= 6; i++) {
      const st = (activeTextConfig.stages && activeTextConfig.stages[i]) || DEFAULT_CANVAS_TEXT_CONFIG.stages[i];
      const isHidden = !!(activeTextConfig._hidden && activeTextConfig._hidden[`stages.${i}`]);
      html += `
        <div class="p-3 rounded-xl bg-slate-950/60 border border-slate-800 space-y-2">
          <div class="flex items-center justify-between">
            <span class="text-xs font-bold text-indigo-400">Stage ${i}</span>
            <label class="flex items-center gap-1.5 cursor-pointer text-xs text-slate-300 select-none">
              <input type="checkbox" id="cfg-stage-${i}-show" ${!isHidden ? 'checked' : ''} class="rounded w-3.5 h-3.5 accent-indigo-500 cursor-pointer">
              <span>Show</span>
            </label>
          </div>
          <div>
            <label class="text-[11px] text-slate-400 block mb-0.5">Title Banner</label>
            <input type="text" id="cfg-stage-${i}-title" value="${escapeHtml(st.title || '')}" class="w-full bg-slate-900 border border-slate-700 rounded-lg px-2.5 py-1 text-xs text-white font-mono outline-none focus:border-indigo-400">
          </div>
          <div>
            <label class="text-[11px] text-slate-400 block mb-0.5">Scientific Description</label>
            <textarea id="cfg-stage-${i}-desc" rows="2" class="w-full bg-slate-900 border border-slate-700 rounded-lg p-2 text-xs text-slate-200 outline-none focus:border-indigo-400 resize-none">${escapeHtml(st.desc || '')}</textarea>
          </div>
        </div>
      `;
    }
    html += '</div>';
    container.innerHTML = html;
    return;
  }

  const fields = CONFIG_FIELD_META[currentCfgTab] || [];
  const sectionData = activeTextConfig[currentCfgTab] || {};
  let html = `
    <div class="flex items-center justify-between pb-2 mb-3 border-b border-slate-800">
      <span class="text-xs text-slate-400 font-medium">Text Labels in Tab</span>
      <div class="flex items-center gap-2">
        <button type="button" onclick="setAllFieldToggles(true)" class="px-2 py-0.5 text-[11px] rounded bg-slate-800 hover:bg-slate-700 text-cyan-300 transition">Show All</button>
        <button type="button" onclick="setAllFieldToggles(false)" class="px-2 py-0.5 text-[11px] rounded bg-slate-800 hover:bg-slate-700 text-slate-400 hover:text-white transition">Hide All</button>
      </div>
    </div>
    <div class="space-y-3">
  `;
  fields.forEach(f => {
    const val = sectionData[f.key] !== undefined ? sectionData[f.key] : "";
    const isHidden = !!(activeTextConfig._hidden && activeTextConfig._hidden[`${currentCfgTab}.${f.key}`]);
    html += `
      <div class="p-2.5 rounded-xl bg-slate-950/60 border border-slate-800 space-y-1.5">
        <div class="flex justify-between items-center">
          <label class="text-xs font-semibold text-slate-200 flex items-center gap-2 cursor-pointer select-none">
            <input type="checkbox" id="cfg-show-${f.key}" ${!isHidden ? 'checked' : ''} class="rounded w-3.5 h-3.5 accent-indigo-500 cursor-pointer">
            <span>${f.label}</span>
          </label>
          ${f.hint ? `<span class="text-[10px] text-cyan-400 font-mono">${f.hint}</span>` : ''}
        </div>
        <input type="text" id="cfg-field-${f.key}" value="${escapeHtml(val)}" class="w-full bg-slate-900 border border-slate-700 rounded-lg px-2.5 py-1.5 text-xs text-white font-mono outline-none focus:border-indigo-400 transition">
      </div>
    `;
  });
  html += '</div>';
  container.innerHTML = html;
}

function setAllFieldToggles(show) {
  const fields = CONFIG_FIELD_META[currentCfgTab] || [];
  fields.forEach(f => {
    const cb = document.getElementById(`cfg-show-${f.key}`);
    if (cb) cb.checked = show;
  });
}
window.setAllFieldToggles = setAllFieldToggles;

function setAllStageToggles(show) {
  for (let i = 1; i <= 6; i++) {
    const cb = document.getElementById(`cfg-stage-${i}-show`);
    if (cb) cb.checked = show;
  }
}
window.setAllStageToggles = setAllStageToggles;

function escapeHtml(str) {
  return String(str).replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/'/g, '&#39;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function applyTextConfigFromModal() {
  if (!activeTextConfig._hidden) activeTextConfig._hidden = {};

  if (currentCfgTab === 'stages') {
    if (!activeTextConfig.stages) activeTextConfig.stages = {};
    for (let i = 1; i <= 6; i++) {
      const titleEl = document.getElementById(`cfg-stage-${i}-title`);
      const descEl = document.getElementById(`cfg-stage-${i}-desc`);
      const showCb = document.getElementById(`cfg-stage-${i}-show`);
      if (titleEl && descEl) {
        activeTextConfig.stages[i] = {
          title: titleEl.value,
          desc: descEl.value
        };
      }
      if (showCb) {
        activeTextConfig._hidden[`stages.${i}`] = !showCb.checked;
      }
    }
  } else {
    const fields = CONFIG_FIELD_META[currentCfgTab] || [];
    if (!activeTextConfig[currentCfgTab]) activeTextConfig[currentCfgTab] = {};
    fields.forEach(f => {
      const el = document.getElementById(`cfg-field-${f.key}`);
      const showCb = document.getElementById(`cfg-show-${f.key}`);
      if (el) {
        activeTextConfig[currentCfgTab][f.key] = el.value;
      }
      if (showCb) {
        activeTextConfig._hidden[`${currentCfgTab}.${f.key}`] = !showCb.checked;
      }
    });
  }

  if (typeof updateStageText === 'function') {
    updateStageText();
  }
  closeTextConfigModal();
}

function resetTextConfigToDefaults() {
  activeTextConfig = deepMerge(DEFAULT_CANVAS_TEXT_CONFIG, (typeof CANVAS_TEXT_CONFIG !== 'undefined') ? CANVAS_TEXT_CONFIG : {});
  activeTextConfig._hidden = {};
  renderConfigForm();
  if (typeof updateStageText === 'function') {
    updateStageText();
  }
}

function exportConfigJs() {
  const content = `/**
 * pyPLP Simulator - Canvas Text & Annotation Configuration
 * =========================================================
 * Generated on: ${new Date().toISOString()}
 */

const CANVAS_TEXT_CONFIG = ${JSON.stringify(activeTextConfig, null, 2)};

if (typeof window !== 'undefined') {
  window.CANVAS_TEXT_CONFIG = CANVAS_TEXT_CONFIG;
}
if (typeof module !== 'undefined' && module.exports) {
  module.exports = CANVAS_TEXT_CONFIG;
}
`;
  const blob = new Blob([content], { type: "application/javascript" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = "canvas_config.js";
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}
