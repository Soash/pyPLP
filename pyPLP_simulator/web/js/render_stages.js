/**
 * pyPLP Simulator - Canvas Reaction Stages Rendering (Stages 1 through 6)
 */

// Vibrant Glowing Biotin reporter tag (with flexible chemical linker like FAM)
function drawBiotinTag(ctx, arg1, arg2, arg3, arg4, arg5) {
  let attachPt = null;
  let x = 0, y = 0, alpha = 1.0, label = "Biotin";

  if (typeof arg1 === 'object' && arg1 !== null) {
    // Signature: (ctx, attachPt, x, y, alpha, label)
    attachPt = arg1;
    x = arg2;
    y = arg3;
    alpha = (arg4 !== undefined) ? arg4 : 1.0;
    label = (arg5 !== undefined) ? arg5 : "Biotin";
  } else {
    // Signature: (ctx, x, y, alpha, label)
    x = arg1;
    y = arg2;
    alpha = (arg3 !== undefined) ? arg3 : 1.0;
    label = (arg4 !== undefined) ? arg4 : "Biotin";
  }

  ctx.save();
  ctx.globalAlpha = alpha;

  // Flexible dashed chemical linker (----) connecting 5' nucleotide to Biotin (like FAM)
  if (attachPt) {
    ctx.save();
    ctx.strokeStyle = "#94a3b8";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([2, 2]);
    ctx.beginPath();
    ctx.moveTo(attachPt.x, attachPt.y);
    ctx.lineTo(x, y);
    ctx.stroke();
    ctx.restore();
  }

  const pulse = Math.sin(performance.now() * 0.008) * 3 + 8; // Glowing pulse like FAM
  
  // Outer pulsing halo
  ctx.fillStyle = "rgba(245, 158, 11, 0.35)";
  ctx.shadowColor = "#f59e0b";
  ctx.shadowBlur = 18 + pulse * 1.5;
  ctx.beginPath();
  ctx.arc(x, y, pulse + 4, 0, Math.PI * 2);
  ctx.fill();

  // Vibrant core sphere with radial gradient
  const grad = ctx.createRadialGradient(x - 2, y - 2, 1, x, y, pulse);
  grad.addColorStop(0, "#fef08a"); // Bright specular highlight
  grad.addColorStop(0.4, "#fbbf24");
  grad.addColorStop(1, "#d97706");
  ctx.fillStyle = grad;
  ctx.strokeStyle = "#b45309";
  ctx.lineWidth = 1.2;
  ctx.beginPath();
  ctx.arc(x, y, pulse, 0, Math.PI * 2);
  ctx.fill();
  ctx.stroke();

  // Label badge (clean "Biotin")
  ctx.shadowBlur = 6;
  ctx.shadowColor = "#000000";
  ctx.fillStyle = "#ffffff";
  ctx.font = "bold 8px 'JetBrains Mono', monospace";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(label, x, y + pulse + 10);

  ctx.restore();
}

// --- STAGE 1: SOLUTION STATE ---
function drawStage1(w, h, p) {
  const cx = w * 0.5;
  const targetY = h * 0.67;
  const spacing = renderSettings.spacing;
  const radius = renderSettings.radius;
  const isBasesVisible = shouldDrawBases();
  const proto = state.protocol || { hasTarget: false, hasPlp: false, hasLigase: false };

  if (!proto.hasTarget && !proto.hasPlp) {
    ctx.save();
    const boxW = Math.min(w * 0.75, 580);
    const boxH = 150;
    const boxX = cx - boxW * 0.5;
    const boxY = h * 0.48 - boxH * 0.5;

    ctx.fillStyle = "rgba(15, 23, 42, 0.75)";
    ctx.strokeStyle = "rgba(56, 189, 248, 0.35)";
    ctx.lineWidth = 1.5;
    ctx.setLineDash([6, 6]);
    ctx.beginPath();
    ctx.roundRect(boxX, boxY, boxW, boxH, 16);
    ctx.fill();
    ctx.stroke();

    ctx.setLineDash([]);
    ctx.fillStyle = "#38bdf8";
    ctx.font = "bold 15px Inter";
    ctx.textAlign = "center";
    ctx.fillText("Empty Reaction Buffer", cx, boxY + 38);

    ctx.fillStyle = "#cbd5e1";
    ctx.font = "12px JetBrains Mono";
    ctx.fillText("Use the Protocol Actions toolbar above to pipette components:", cx, boxY + 66);

    ctx.fillStyle = "#94a3b8";
    ctx.font = "11px JetBrains Mono";
    ctx.fillText("1. Click '+ Add Target' to load target sequence into solution", cx, boxY + 90);
    ctx.fillText("2. Click '+ Add PLP' to load padlock probe into solution", cx, boxY + 110);
    ctx.fillText("3. Click '+ Add Ligase' to initiate covalent ligation", cx, boxY + 130);
    ctx.restore();
    return;
  }

  // 1. TARGET STRAND (if added to solution)
  if (proto.hasTarget) {
    const targetLen = (state.target && state.target.length) || 50;
    const targetW = (targetLen - 1) * spacing;
    const startX = cx - targetW * 0.5;
    
    const targetPts = [];
    for (let i = 0; i < targetLen; i++) {
      const x = startX + i * spacing;
      const y = targetY + Math.sin(i * 0.3 + p * 2) * 4;
      targetPts.push({ x, y });
    }
    drawStrandWithSequence(targetPts, state.target, "#f59e0b", getCanvasText('strands', 'targetStrandLabel', "Target ({targetType}) 5' → 3'", { targetType: state.targetType }), false);
  }

  // 2. FULL LINEAR PADLOCK PROBE (SINGLE CONTINUOUS STRAND)
  const plpY = h * 0.32 + Math.sin(p * 3) * 8;
  const lLen = (state.leftArmSeq && state.leftArmSeq.length) || 15;
  const rLen = (state.rightArmSeq && state.rightArmSeq.length) || 15;
  const bbLen = (state.backboneSeq && state.backboneSeq.length) || 40;
  const fullPlpSeq = state.plp || (state.leftArmSeq + state.backboneSeq + state.rightArmSeq);
  const totalPlpLen = fullPlpSeq.length;

  if (proto.hasPlp) {
    const totalPlpW = (totalPlpLen - 1) * spacing;
    const plpStartX = cx - totalPlpW * 0.5;

    const plpPts = [];
    for (let i = 0; i < totalPlpLen; i++) {
      plpPts.push({
        x: plpStartX + i * spacing,
        y: plpY + Math.sin(i * 0.18 + p * 2.5) * 6
      });
    }

    // Draw Backbone Line across the entire continuous strand
    ctx.save();
    ctx.strokeStyle = isBasesVisible ? "rgba(148, 163, 184, 0.35)" : "#a855f7";
    ctx.lineWidth = isBasesVisible ? 2.5 : 4;
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.beginPath();
    ctx.moveTo(plpPts[0].x, plpPts[0].y);
    for (let i = 1; i < plpPts.length; i++) {
      ctx.lineTo(plpPts[i].x, plpPts[i].y);
    }
    ctx.stroke();

    // Draw domain-colored segments when zoomed out
    if (!isBasesVisible) {
      // Left Arm segment
      ctx.strokeStyle = "#06b6d4";
      ctx.lineWidth = 4;
      ctx.beginPath();
      ctx.moveTo(plpPts[0].x, plpPts[0].y);
      for (let i = 1; i < lLen; i++) ctx.lineTo(plpPts[i].x, plpPts[i].y);
      ctx.stroke();

      // Right Arm segment
      ctx.strokeStyle = "#10b981";
      ctx.lineWidth = 4;
      ctx.beginPath();
      ctx.moveTo(plpPts[totalPlpLen - rLen].x, plpPts[totalPlpLen - rLen].y);
      for (let i = totalPlpLen - rLen + 1; i < totalPlpLen; i++) ctx.lineTo(plpPts[i].x, plpPts[i].y);
      ctx.stroke();

      // Ticks along strand
      for (let i = 0; i < plpPts.length; i += 2) {
        ctx.fillStyle = i < lLen ? "#06b6d4" : (i >= totalPlpLen - rLen ? "#10b981" : "#a855f7");
        ctx.beginPath();
        ctx.arc(plpPts[i].x, plpPts[i].y, 2.5, 0, Math.PI * 2);
        ctx.fill();
      }
    } else {
      // When zoomed in: RENDER EVERY SINGLE NUCLEOTIDE OF FULL PLP (Left Arm, Backbone, Right Arm)
      for (let i = 0; i < totalPlpLen; i++) {
        const pt = plpPts[i];
        const base = fullPlpSeq[i];
        let sub = '';
        if (i === 0) sub = "5'PO4";
        else if (i === lLen - 1) sub = "Left";
        else if (i === totalPlpLen - rLen) sub = "Right";
        else if (i === totalPlpLen - 1) sub = "3'OH";
        drawNucleotideNode(ctx, pt.x, pt.y, base, radius, false, sub);
      }
    }

    // Domain Banner Annotations above the full probe
    ctx.save();
    ctx.font = "bold 11px JetBrains Mono";
    
    // Left Arm label
    ctx.fillStyle = "#06b6d4";
    const leftMid = plpPts[Math.floor(lLen * 0.5)];
    ctx.textAlign = "center";
    ctx.fillText(getCanvasText('strands', 'leftArmLinearLabel', "5'-Left Arm"), leftMid.x, leftMid.y - 20);

    // Backbone label
    ctx.fillStyle = "#c084fc";
    const bbMid = plpPts[Math.floor(lLen + bbLen * 0.5)];
    ctx.fillText(getCanvasText('strands', 'backboneLinearLabel', "Backbone"), bbMid.x, bbMid.y - 20);

    // Right Arm label
    ctx.fillStyle = "#10b981";
    const rightMid = plpPts[Math.floor(totalPlpLen - rLen * 0.5)];
    ctx.fillText(getCanvasText('strands', 'rightArmLinearLabel', "Right Arm-3'"), rightMid.x, rightMid.y - 20);

    // Overall Continuous Strand Title
    ctx.fillStyle = "#e2e8f0";
    ctx.font = "bold 11px JetBrains Mono";
    ctx.fillText("Full Linear Padlock Probe (PLP 5' → 3')", cx, plpY - 42);
    ctx.restore();
    ctx.restore();
  } else {
    // Subtle placeholder when PLP is not loaded
    ctx.save();
    ctx.strokeStyle = "rgba(168, 85, 247, 0.2)";
    ctx.setLineDash([4, 6]);
    ctx.lineWidth = 1.5;
    ctx.strokeRect(cx - 200, plpY - 20, 400, 40);
    ctx.fillStyle = "rgba(168, 85, 247, 0.4)";
    ctx.font = "italic 11px JetBrains Mono";
    ctx.textAlign = "center";
    ctx.fillText("PLP empty: Click '+ Add PLP' to load Padlock Probe", cx, plpY + 4);
    ctx.restore();
  }

  // 3. LIGASE ENZYME (if added to solution)
  if (proto.hasLigase) {
    const ligX = cx + Math.min(w * 0.35, 240);
    const ligY = plpY - 35 + Math.sin(p * 2 + 1) * 7;
    ctx.save();
    ctx.fillStyle = "rgba(16, 185, 129, 0.2)";
    ctx.beginPath();
    ctx.arc(ligX, ligY, 32, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#10b981";
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 12;
    ctx.beginPath();
    ctx.arc(ligX, ligY, 18, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 10px JetBrains Mono";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    const enzymeName = state.targetType === 'RNA' 
      ? getCanvasText('ligation', 'enzymeRnaLigase', "SplintR") 
      : getCanvasText('ligation', 'enzymeDnaLigase', "Ligase");
    ctx.fillText(enzymeName, ligX, ligY);

    ctx.font = "bold 9px JetBrains Mono";
    ctx.fillStyle = "#6ee7b7";
    ctx.fillText("Enzyme Ready", ligX, ligY + 28);
    ctx.restore();
  }
}

// --- STAGE 2: TARGET HYBRIDIZATION ---
function drawStage2(w, h, p) {
  const geo = computeHorseshoeGeometry(w, h, p);
  const isBasesVisible = geo.isBasesVisible;
  const radius = geo.radius;

  // Draw Target Strand
  drawStrandWithSequence(geo.targetPts, state.target, "#f59e0b", getCanvasText('strands', 'targetStrandLabel', "Target ({targetType}) 5' → 3'", { targetType: state.targetType }), false);

  // Draw Left Arm & Right Arm (Left label at outer left tip, Right label at outer right tip)
  drawStrandWithSequence(geo.lArmPts, state.leftArmSeq, "#06b6d4", getCanvasText('strands', 'leftArmHybridLabel', "Left Arm (5'-PO4)"), true);
  drawStrandWithSequence(geo.rArmPts, state.rightArmSeq, "#10b981", getCanvasText('strands', 'rightArmHybridLabel', "Right Arm (3'-OH)"), false);

  // Backbone Arched Loop
  ctx.save();
  ctx.strokeStyle = isBasesVisible ? "rgba(168, 85, 247, 0.45)" : "#a855f7";
  ctx.lineWidth = isBasesVisible ? 2.5 : 4;
  ctx.shadowColor = "#a855f7";
  ctx.shadowBlur = isBasesVisible ? 0 : 10;
  ctx.beginPath();
  ctx.moveTo(geo.bbStart.x, geo.bbStart.y);
  ctx.bezierCurveTo(geo.cp1.x, geo.cp1.y, geo.cp2.x, geo.cp2.y, geo.bbEnd.x, geo.bbEnd.y);
  ctx.stroke();

  if (isBasesVisible && state.backboneSeq) {
    for (let i = 0; i < geo.bbLen; i++) {
      const pt = geo.bbPts[i];
      const base = state.backboneSeq[i];
      drawNucleotideNode(ctx, pt.x, pt.y, base, radius, false, '');
    }
  } else {
    ctx.fillStyle = "#c084fc";
    for (let i = 0; i < geo.bbPts.length; i += 3) {
      ctx.beginPath();
      ctx.arc(geo.bbPts[i].x, geo.bbPts[i].y, 2.5, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  if (geo.loopHeight > 25) {
    ctx.fillStyle = "#c084fc";
    ctx.font = "bold 11px JetBrains Mono";
    ctx.textAlign = "center";
    ctx.fillText(getCanvasText('strands', 'backboneLoopLabel', "Backbone Loop"), (geo.bbStart.x + geo.bbEnd.x) * 0.5, (geo.bbStart.y + geo.bbEnd.y) * 0.5 - geo.loopHeight * 0.82);
  }
  ctx.restore();

  // Hydrogen Bond Rungs between Arm and Target (appear as arms land on target)
  if (p > 0.60) {
    ctx.save();
    const topY = geo.currentArmY + (isBasesVisible ? radius + 2 : 4);
    const botY = geo.targetY - (isBasesVisible ? radius + 2 : 4);

    // Left Arm Base Pairing
    for (let k = 0; k < geo.lLen; k++) {
      const pt = geo.lArmPts[k];
      const armBase = state.leftArmSeq[k];
      let targetBase = 'N';
      if (state.leftTargetIdx !== -1) {
        const tIdx = state.leftTargetIdx + (geo.lLen - 1 - k);
        if (tIdx >= 0 && tIdx < state.target.length) targetBase = state.target[tIdx];
      }
      const matched = isWatsonCrick(armBase, targetBase);
      ctx.strokeStyle = matched ? "rgba(52, 211, 153, 0.85)" : "rgba(244, 63, 94, 0.9)";
      ctx.lineWidth = matched ? 2 : 2.5;
      ctx.setLineDash(matched ? [2, 3] : []);
      ctx.beginPath();
      ctx.moveTo(pt.x, topY);
      ctx.lineTo(pt.x, botY);
      ctx.stroke();
    }

    // Right Arm Base Pairing
    for (let j = 0; j < geo.rLen; j++) {
      const pt = geo.rArmPts[j];
      const armBase = state.rightArmSeq[j];
      let targetBase = 'N';
      if (state.rightTargetIdx !== -1) {
        const tIdx = state.rightTargetIdx + (geo.rLen - 1 - j);
        if (tIdx >= 0 && tIdx < state.target.length) targetBase = state.target[tIdx];
      }
      const matched = isWatsonCrick(armBase, targetBase);
      ctx.strokeStyle = matched ? "rgba(52, 211, 153, 0.85)" : "rgba(244, 63, 94, 0.9)";
      ctx.lineWidth = matched ? 2 : 2.5;
      ctx.setLineDash(matched ? [2, 3] : []);
      ctx.beginPath();
      ctx.moveTo(pt.x, topY);
      ctx.lineTo(pt.x, botY);
      ctx.stroke();
    }
    ctx.restore();
  }

  // Ligation Nick Callout / Proofreading Indicator
  if (p > 0.7) {
    ctx.save();
    const calloutColor = state.isLigatable ? "#10b981" : "#f43f5e";
    ctx.strokeStyle = calloutColor;
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(geo.nickX + geo.gapOffset * 0.5, geo.landingY, 15, 0, Math.PI * 2);
    ctx.stroke();

    ctx.fillStyle = calloutColor;
    ctx.font = "bold 11px JetBrains Mono";
    ctx.textAlign = "center";
    let msg = getCanvasText('ligation', 'nickFlushBadge', "0-bp Nick (3'-OH | 5'-PO4)");
    if (!state.isLigatable) {
      if (state.hasNickMismatch) {
        msg = getCanvasText('ligation', 'nickMismatchBadge', "3'-End Mismatch (Arrested)");
      } else if (state.nickGap !== 999) {
        msg = getCanvasText('ligation', 'nickGapBadge', "{gap}-bp Gap (Arrested)", { gap: state.nickGap });
      } else {
        msg = getCanvasText('ligation', 'nickNoAlignBadge', "No Alignment");
      }
    }
    ctx.fillText(msg, geo.nickX + geo.gapOffset * 0.5, geo.landingY - (isBasesVisible ? radius + 18 : 24));
    ctx.restore();
  }

  // Ligase presence in buffer during Stage 2 hybridization
  if (state.currentStage === 2 && state.protocol && state.protocol.hasLigase) {
    ctx.save();
    const ligX = geo.nickX + 140;
    const ligY = geo.landingY - 80;
    ctx.fillStyle = "rgba(16, 185, 129, 0.2)";
    ctx.beginPath();
    ctx.arc(ligX, ligY, 24, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#10b981";
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 12;
    ctx.beginPath();
    ctx.arc(ligX, ligY, 15, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 9px 'JetBrains Mono', monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    const enzymeName = state.targetType === 'RNA' ? "SplintR" : "Ligase";
    ctx.fillText(enzymeName, ligX, ligY);

    ctx.font = "bold 8px 'JetBrains Mono', monospace";
    ctx.fillStyle = "#6ee7b7";
    ctx.fillText("Ready in Buffer", ligX, ligY + 22);
    ctx.restore();
  }
}

// --- STAGE 3: LIGATION & TOPOLOGY ---
function drawStage3(w, h, p) {
  drawStage2(w, h, 1.0);

  const geo = computeHorseshoeGeometry(w, h, 1.0);
  const ligBufferX = geo.nickX + 140;
  const ligBufferY = geo.landingY - 80;

  // Smooth seamless glide from buffer position (end of Stage 2) to the nick during p in [0.0, 0.45]
  const t = Math.min(p / 0.45, 1.0);
  const ease = t * t * (3 - 2 * t); // smoothstep
  const enzymeX = ligBufferX + (geo.nickX - ligBufferX) * ease;
  const enzymeY = ligBufferY + (geo.landingY - ligBufferY) * ease;
  ctx.save();

  if (state.isLigatable) {
    ctx.fillStyle = "rgba(16, 185, 129, 0.25)";
    ctx.beginPath();
    ctx.arc(enzymeX, enzymeY, 28, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#10b981";
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 15;
    ctx.beginPath();
    ctx.arc(enzymeX, enzymeY, 15, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 9px Inter";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(state.targetType === 'RNA' ? getCanvasText('ligation', 'enzymeRnaLigase', "SplintR") : getCanvasText('ligation', 'enzymeDnaLigase', "Ligase"), enzymeX, enzymeY);

    // Status annotation above/below enzyme
    ctx.font = "bold 8px 'JetBrains Mono', monospace";
    if (p < 0.45) {
      ctx.fillStyle = "#6ee7b7";
      ctx.fillText("Approaching Nick...", enzymeX, enzymeY + 22);
    } else if (p < 0.6) {
      ctx.fillStyle = "#34d399";
      ctx.fillText("Interrogating Nick...", enzymeX, enzymeY - 22);
    }

    if (p > 0.5) {
      ctx.strokeStyle = "#34d399";
      ctx.lineWidth = 5;
      ctx.shadowColor = "#34d399";
      ctx.shadowBlur = 25;
      ctx.beginPath();
      ctx.moveTo(geo.nickX - geo.spacing * 0.5, geo.landingY);
      ctx.lineTo(geo.nickX + geo.spacing * 0.5, geo.landingY);
      ctx.stroke();

      ctx.fillStyle = "#34d399";
      ctx.font = "bold 12px JetBrains Mono";
      ctx.textAlign = "center";
      ctx.fillText(getCanvasText('ligation', 'sealedSuccessBanner', "COVALENT PHOSPHODIESTER BOND SEALED"), geo.cx, geo.landingY + 45);
    }
  } else {
    ctx.fillStyle = "rgba(244, 63, 94, 0.25)";
    ctx.beginPath();
    ctx.arc(enzymeX, enzymeY, 28, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#f43f5e";
    ctx.shadowColor = "#f43f5e";
    ctx.shadowBlur = 15;
    ctx.beginPath();
    ctx.arc(enzymeX, enzymeY, 15, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 13px Inter";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText("✕", enzymeX, enzymeY);

    ctx.fillStyle = "#f43f5e";
    ctx.font = "bold 12px JetBrains Mono";
    ctx.textAlign = "center";
    const errMsg = state.hasNickMismatch 
      ? getCanvasText('ligation', 'arrestedMismatchBanner', "LIGATION ARRESTED: 3' MISMATCH DETECTED") 
      : (state.nickGap !== 999 ? getCanvasText('ligation', 'arrestedGapBanner', "LIGATION ARRESTED: {gap}-BP GAP DETECTED", { gap: state.nickGap }) : getCanvasText('ligation', 'arrestedNoAlignBanner', "LIGATION ARRESTED: ARMS DO NOT ALIGN"));
    ctx.fillText(errMsg, geo.cx, geo.landingY + 45);
  }
  ctx.restore();
}

// --- STAGE 4: RCA PRIMER ANNEALING ON HYBRIDIZED LOOP (SAME SCENE) ---
function drawStage4(w, h, p) {
  // 1. Draw hybridized padlock on target with sealed nick (identical to Stage 3!)
  drawStage3(w, h, 1.0);

  const geo = computeHorseshoeGeometry(w, h, 1.0);
  const isBasesVisible = geo.isBasesVisible;
  const radius = geo.radius;
  const primerLen = (state.primer && state.primer.length) || 18;
  const total = geo.totalPlpLen;

  if (state.protocol && state.protocol.hasPrimer === false) {
    ctx.save();
    ctx.fillStyle = "rgba(244, 63, 94, 0.4)";
    ctx.font = "italic 11px JetBrains Mono";
    ctx.textAlign = "center";
    ctx.fillText("Primer empty: Click '+ Add RCA Primer' to load primer", geo.cx, geo.landingY - geo.loopHeight - 35);
    ctx.restore();
    return;
  }

  if (state.isPrimerComplementary) {
    // COMPLEMENTARY PRIMER ANNEALS ANTIPARALLEL DIRECTLY ONTO THE BACKBONE LOOP!
    let dockStartIdx = state.primerDockIdx;
    if (dockStartIdx === -1 || dockStartIdx < 0 || dockStartIdx >= total) {
      dockStartIdx = Math.floor(geo.lLen + geo.bbLen * 0.28);
    }

    const isDocked = (state.protocol && state.protocol.isPrimerAnnealed) || p > 0.6;
    const dockDist = isDocked ? 0 : (1 - Math.min(p * 1.3, 1.0)) * 60;

    // Primer is 5' -> 3' (k = 0 is 5', k = primerLen - 1 is 3'-OH).
    // Template runs 5' -> 3' (left to right).
    // Because DNA is strictly ANTIPARALLEL:
    // Primer base k (5' -> 3') pairs with template at (dockStartIdx + (primerLen - 1 - k))!
    // - Primer 5' (k = 0) is on the RIGHT.
    // - Primer 3'-OH (k = primerLen - 1) is on the LEFT.
    const pPts = [];
    const templatePts = [];

    for (let k = 0; k < primerLen; k++) {
      const plpIdx = (dockStartIdx + (primerLen - 1 - k)) % total;
      const bPt = geo.plpPts[plpIdx] || geo.plpPts[geo.lLen];
      templatePts.push(bPt);

      // Normal outward from loop apex center
      const dx = bPt.x - geo.cx;
      const dy = bPt.y - (geo.landingY - geo.loopHeight * 0.45);
      const len = Math.hypot(dx, dy) || 1;
      const nx = dx / len;
      const ny = dy / len;

      const px = bPt.x + nx * (dockDist + (isBasesVisible ? radius * 2 + 14 : 16));
      const py = bPt.y + ny * (dockDist + (isBasesVisible ? radius * 2 + 14 : 16));
      pPts.push({ x: px, y: py });
    }

    drawStrandWithSequence(pPts, state.primer, "#f43f5e", getCanvasText('rcaPrimer', 'primerStrandLabel', "RCA Primer (3'-OH ◄── 5')"), false);

    // Terminal Badges: 5' on the right, 3'-OH on the left
    if (pPts.length > 0) {
      ctx.save();
      ctx.font = "bold 9px JetBrains Mono";
      // 5' terminal badge at pPts[0] (right)
      ctx.fillStyle = "#f43f5e";
      ctx.textAlign = "center";
      ctx.fillText("5'", pPts[0].x + 14, pPts[0].y - 8);

      // Glowing Biotin reporter tag tethered via flexible chemical linker (----)
      drawBiotinTag(ctx, pPts[0], pPts[0].x + 22, pPts[0].y - 20, 1.0, "Biotin");

      // 3'-OH terminal badge at pPts[primerLen - 1] (left)
      ctx.fillStyle = "#38bdf8";
      ctx.font = "bold 9px JetBrains Mono";
      ctx.fillText("3'-OH", pPts[primerLen - 1].x - 18, pPts[primerLen - 1].y - 8);
      ctx.restore();
    }

    // Base-pairing rungs when docked onto the loop (strictly Watson-Crick verified)
    if (p > 0.4 || isDocked) {
      ctx.save();
      for (let k = 0; k < primerLen; k++) {
        const plpIdx = (dockStartIdx + (primerLen - 1 - k)) % total;
        const pBase = state.primer[k];
        const tBase = (state.plp && state.plp[plpIdx]) || 'N';
        const matched = isWatsonCrick(pBase, tBase);

        ctx.strokeStyle = matched ? "rgba(52, 211, 153, 0.9)" : "rgba(244, 63, 94, 0.9)";
        ctx.lineWidth = matched ? 1.8 : 2.0;
        ctx.setLineDash(matched ? [2, 3] : []);
        ctx.beginPath();
        ctx.moveTo(pPts[k].x, pPts[k].y);
        ctx.lineTo(templatePts[k].x, templatePts[k].y);
        ctx.stroke();
      }
      ctx.restore();

      // Annealed Banner
      ctx.save();
      ctx.fillStyle = "#34d399";
      ctx.font = "bold 12px JetBrains Mono";
      ctx.textAlign = "center";
      ctx.fillText(getCanvasText('rcaPrimer', 'primerAnnealedBanner', "RCA PRIMER ANNEALED TO BACKBONE (3'-OH EXTENSION READY)"), geo.cx, geo.landingY - geo.loopHeight - 48);
      ctx.restore();
    }
  } else {
    // NON-COMPLEMENTARY PRIMER: floats unannealed above loop
    const floatDist = 70 + Math.sin(p * 3.5) * 8;
    const pPts = [];
    const startX = geo.cx - ((primerLen - 1) * geo.spacing) * 0.5;
    const floatY = geo.landingY - geo.loopHeight - floatDist;
    for (let k = 0; k < primerLen; k++) {
      pPts.push({ x: startX + k * geo.spacing, y: floatY + Math.sin(k * 0.3) * 4 });
    }

    drawStrandWithSequence(pPts, state.primer, "#f43f5e", getCanvasText('rcaPrimer', 'primerStrandLabel', "5'-RCA Primer (Unbound)-3'OH"), false);

    // Warning Badge
    ctx.save();
    const midPt = pPts[Math.floor(primerLen * 0.5)];
    ctx.fillStyle = "rgba(244, 63, 94, 0.25)";
    ctx.strokeStyle = "#f43f5e";
    ctx.lineWidth = 1.5;
    const badgeW = 210;
    ctx.beginPath();
    ctx.roundRect(midPt.x - badgeW * 0.5, midPt.y - 34, badgeW, 22, [6]);
    ctx.fill();
    ctx.stroke();

    ctx.fillStyle = "#ffffff";
    ctx.font = "bold 10px JetBrains Mono";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(getCanvasText('rcaPrimer', 'primerBlockedBadge', "✕ NON-COMPLEMENTARY PRIMER"), midPt.x, midPt.y - 23);

    ctx.fillStyle = "#f43f5e";
    ctx.font = "bold 12px JetBrains Mono";
    ctx.fillText(getCanvasText('rcaPrimer', 'primerArrestedHeader', "PRIMER ANNEALING ARRESTED: NO COMPLEMENTARY BINDING SITE"), geo.cx, geo.landingY + 68);
    ctx.font = "10px JetBrains Mono";
    ctx.fillStyle = "#fda4af";
    ctx.fillText(getCanvasText('rcaPrimer', 'primerArrestedDetail', "Primer sequence matches 0 cassettes on circular template ({matchPct}% match < 70% required)", { matchPct: state.primerMatchPct }), geo.cx, geo.landingY + 86);
    ctx.restore();
  }
}

// --- STAGE 5: PHI29 ENGAGEMENT, TARGET DISPLACEMENT & CONTINUOUS RCA ---
function drawStage5(w, h, p, dt) {
  const canRcaRun = state.isLigatable && state.isPrimerComplementary;
  const spacing = renderSettings.spacing;
  const isBasesVisible = shouldDrawBases();
  const userRepeats = Math.max(1, (state.protocol && state.protocol.userRepeats) || 3);

  if (!canRcaRun) {
    drawStage4(w, h, 1.0);
    state.concatemerRepeats = 0;
    const repCounter = document.getElementById('repeats-counter');
    if (repCounter) repCounter.innerText = "0";

    ctx.save();
    ctx.fillStyle = "#f43f5e";
    ctx.font = "bold 13px JetBrains Mono";
    ctx.textAlign = "center";
    const reasonKey = !state.isLigatable ? 'rcaReasonLigation' : 'rcaReasonPrimer';
    const reasonDefault = !state.isLigatable 
      ? "LIGATION ARRESTED (PADLOCK NOT CIRCULARIZED)" 
      : "RCA PRIMER NOT ANNEALED (NO 3'-OH TERMINUS)";
    const reason = getCanvasText('rollingCircle', reasonKey, reasonDefault);
    ctx.fillText(getCanvasText('rollingCircle', 'rcaBlockedHeader', `RCA INITIATION BLOCKED: ${reason}`, { reason }), w * 0.5, h * 0.85);
    ctx.font = "10px JetBrains Mono";
    ctx.fillStyle = "#fda4af";
    const subMsg = !state.isLigatable
      ? "Phi29 rolls off linear DNA; covalent circularization is mandatory for rolling circle amplification."
      : "Phi29 DNA Polymerase requires an annealed 3'-OH primer terminus to initiate strand synthesis.";
    ctx.fillText(subMsg, w * 0.5, h * 0.85 + 20);
    ctx.restore();
    return;
  }

  // Compute geometries
  const geoHs = computeHorseshoeGeometry(w, h, 1.0);
  const total = geoHs.totalPlpLen;
  const lLen = geoHs.lLen;
  const bbLen = geoHs.bbLen;
  const rLen = geoHs.rLen;
  const primerLen = (state.primer && state.primer.length) || 18;

  let dockStartIdx = state.primerDockIdx;
  if (dockStartIdx === -1 || dockStartIdx < 0 || dockStartIdx >= total) {
    dockStartIdx = Math.floor(lLen + bbLen * 0.28);
  }
  const dockEndIdx = (dockStartIdx + primerLen - 1) % total;

  // --- REPLICATION STATION CONFIGURATION ---
  // Station on upper-right shoulder of circular template (~1:30 - 2:00 clock position)
  const stationAngle = -Math.PI * 0.22;
  const stationTargetFrac = ((stationAngle - Math.PI * 0.5 + Math.PI * 4) % (Math.PI * 2)) / (Math.PI * 2);
  const stationIdx = Math.round(stationTargetFrac * total) % total;
  const stepsToStation = (dockStartIdx - stationIdx + total) % total;

  // Primer 5' (with Biotin) starting angle when rotAngle = 0:
  const primer5Idx = (dockStartIdx + primerLen - 1) % total;
  const p5BaseAngle = (primer5Idx / total) * Math.PI * 2 + Math.PI * 0.5;
  let rotToStation = (stationAngle - p5BaseAngle) % (Math.PI * 2);
  if (rotToStation < 0) rotToStation += Math.PI * 2;
  if (rotToStation < 0.05) rotToStation += Math.PI * 2;

  // Rotation mechanics:
  // p in [0.0, 0.42]: rotAngle = 0 (Horseshoe & Phi29 displacing target strand)
  // p in [0.42, 0.56]: Target displaced, template morphs smoothly to circle
  // p >= 0.56: Continuous isothermal rolling replication with Phi29 stationary (1:1 with concatemer synthesis)
  // p >= 1.0: RCA repeats completed, template rotation stops
  const isRcaComplete = p >= 1.0;
  let effectiveRotAngle = 0;
  if (p < 0.42) {
    state.rcaRotationAngle = 0;
    effectiveRotAngle = 0;
  } else if (p < 0.56) {
    state.rcaRotationAngle = 0;
    const rotProgress = (p - 0.42) / 0.14;
    effectiveRotAngle = rotProgress * rotToStation;
  } else {
    // Phase 4: Circular template rotation angle strictly matches concatemer synthesis 1:1
    const rcaProgress = Math.min((p - 0.56) / 0.44, 1.0);
    const totalRotation = userRepeats * Math.PI * 2;
    state.rcaRotationAngle = rcaProgress * totalRotation;
    effectiveRotAngle = rotToStation + state.rcaRotationAngle;
  }

  const geoCirc = computeCircularGeometry(w, h, effectiveRotAngle);

  // --- TIMELINE PHASES ---
  // Phase 1 (p in [0.0, 0.10]): Docking onto primer 3'-OH
  // Phase 2 (p in [0.10, 0.42]): Traveling along template & Displacing Target strand base-by-base
  // Phase 3 (p in [0.42, 0.56]): Target displaced & floating away • Circular Template morphing smoothly to station
  // Phase 4 (p >= 0.56): Stationary Phi29 at station • Circular template rotating • Concatemer dynamic spooling

  // Overall uniform morph progress from Horseshoe (0) to Circle (1)
  // CRITICAL: Uniform scalar for ALL points to prevent any backbone tearing or zigzag glitches!
  let m = 0;
  if (p < 0.42) {
    m = 0;
  } else if (p < 0.56) {
    const rawM = (p - 0.42) / 0.14;
    m = rawM * rawM * (3 - 2 * rawM); // smoothstep 0 to 1
  } else {
    m = 1.0;
  }

  // Displacement steps calculation:
  // Primer 3'-OH is at dockStartIdx. Phi29 moves backward: dockStartIdx -> 0 (Left Arm) -> total-1 -> total-rLen (Right Arm)
  const totalDisplaceSteps = dockStartIdx + 1 + rLen;
  const totalStepsToPrimer5 = total - primerLen;
  let stepsTraveled = 0;
  if (p < 0.10) {
    stepsTraveled = 0;
  } else if (p < 0.42) {
    const u = (p - 0.10) / 0.32;
    stepsTraveled = u * totalDisplaceSteps;
  } else if (p < 0.56) {
    const u = (p - 0.42) / 0.14;
    stepsTraveled = totalDisplaceSteps + u * (totalStepsToPrimer5 - totalDisplaceSteps);
  } else {
    stepsTraveled = totalStepsToPrimer5;
  }

  // 1. MORPHED PADLOCK PROBE COORDINATES (UNIFORM CONTINUOUS TOPOLOGY)
  const plpCurrentPts = [];
  for (let i = 0; i < total; i++) {
    const p1 = geoHs.plpPts[i];
    const p2 = geoCirc.plpPts[i];
    plpCurrentPts.push({
      x: (1 - m) * p1.x + m * p2.x,
      y: (1 - m) * p1.y + m * p2.y
    });
  }

  // Draw continuous padlock probe
  ctx.save();
  ctx.strokeStyle = isBasesVisible ? "rgba(168, 85, 247, 0.45)" : "#a855f7";
  ctx.lineWidth = isBasesVisible ? 2.5 : 5;
  ctx.shadowColor = "#a855f7";
  ctx.shadowBlur = isBasesVisible ? 0 : 10;
  ctx.beginPath();
  ctx.moveTo(plpCurrentPts[0].x, plpCurrentPts[0].y);
  for (let i = 1; i < total; i++) {
    ctx.lineTo(plpCurrentPts[i].x, plpCurrentPts[i].y);
  }
  ctx.closePath();
  ctx.stroke();

  const currentRadius = (1 - m) * geoHs.radius + m * geoCirc.bRadius;
  const baseDist = isBasesVisible ? currentRadius * 2 + 14 : 16;
  if (isBasesVisible && state.plp) {
    for (let i = 0; i < total; i++) {
      drawNucleotideNode(ctx, plpCurrentPts[i].x, plpCurrentPts[i].y, state.plp[i], currentRadius);
    }
  } else {
    ctx.fillStyle = "#c084fc";
    for (let i = 0; i < total; i += 2) {
      ctx.beginPath();
      ctx.arc(plpCurrentPts[i].x, plpCurrentPts[i].y, 2.5, 0, Math.PI * 2);
      ctx.fill();
    }
  }
  ctx.restore();

  // Top title
  ctx.save();
  ctx.fillStyle = "#e2e8f0";
  ctx.font = "bold 11px JetBrains Mono";
  ctx.textAlign = "center";
  const topPtY = (1 - m) * (geoHs.landingY - geoHs.loopHeight) + m * (geoCirc.cy - geoCirc.radius);
  let stageTitleText = "Padlock Probe Template";
  if (p < 0.10) stageTitleText = "Phi29 Docking onto Primer 3'-OH";
  else if (p < 0.42) stageTitleText = "⚡ Phi29 Displacing Target Strand Base-by-Base";
  else if (p < 0.56) stageTitleText = "⚡ Target Displaced • Template Morphing to Replication Station";
  else stageTitleText = getCanvasText('rollingCircle', 'stage5Title', "Continuous Rolling Circle Amplification (RCA)");
  const titleOffsetY = (p >= 0.42) ? -50 - baseDist : -45;
  ctx.fillText(stageTitleText, (1 - m) * geoHs.cx + m * geoCirc.cx, topPtY + titleOffsetY);
  ctx.restore();

  // 3. RCA PRIMER COORDINATES ON TEMPLATE (ANTIPARALLEL WITH OUTWARD OFFSET)
  const currentCenter = {
    x: (1 - m) * geoHs.cx + m * geoCirc.cx,
    y: (1 - m) * (geoHs.landingY - geoHs.loopHeight * 0.45) + m * geoCirc.cy
  };

  const primerPts = [];
  const templatePts = [];

  for (let k = 0; k < primerLen; k++) {
    const tIdx = (dockStartIdx + (primerLen - 1 - k)) % total;
    const templatePt = plpCurrentPts[tIdx];
    templatePts.push(templatePt);

    // Vector pointing outward from template loop/circle center
    const dx = templatePt.x - currentCenter.x;
    const dy = templatePt.y - currentCenter.y;
    const len = Math.hypot(dx, dy) || 1;
    const nx = dx / len;
    const ny = dy / len;

    const px = templatePt.x + nx * baseDist;
    const py = templatePt.y + ny * baseDist;
    primerPts.push({ x: px, y: py });
  }

  // 4. PHI29 ENZYME POSITIONING
  const stationPos = {
    x: geoCirc.cx + (geoCirc.radius + baseDist * 0.45) * Math.cos(stationAngle),
    y: geoCirc.cy + (geoCirc.radius + baseDist * 0.45) * Math.sin(stationAngle)
  };

  function getTemplateOuterPos(tIdx, fracNext = 0) {
    const nextIdx = (tIdx - 1 + total) % total;
    const pA = plpCurrentPts[tIdx];
    const pB = plpCurrentPts[nextIdx];
    const tX = (1 - fracNext) * pA.x + fracNext * pB.x;
    const tY = (1 - fracNext) * pA.y + fracNext * pB.y;

    if (m < 0.5 && tY >= geoHs.landingY - 14) {
      return { x: tX, y: tY + (baseDist * 0.45) };
    } else {
      const dx = tX - currentCenter.x;
      const dy = tY - currentCenter.y;
      const len = Math.hypot(dx, dy) || 1;
      return {
        x: tX + (dx / len) * (baseDist * 0.45),
        y: tY + (dy / len) * (baseDist * 0.45)
      };
    }
  }

  let phi29Pos = { x: 0, y: 0 };
  if (p < 0.10) {
    // --- PHASE 1: DOCKING ON PRIMER 3'-OH ---
    const u = p / 0.10;
    const primer3Pt = primerPts[primerLen - 1]; // Primed 3'-OH node of primer!
    phi29Pos = {
      x: primer3Pt.x,
      y: primer3Pt.y - (1 - u) * 60
    };
  } else if (p < 0.56) {
    // --- PHASE 2 & 3: CONTINUOUS TRAVERSAL ALONG TEMPLATE (NO DISCONNECT OR JUMP) ---
    const currentIdx = (dockStartIdx - Math.floor(stepsTraveled) + total * 10) % total;
    const fracStep = stepsTraveled - Math.floor(stepsTraveled);
    phi29Pos = getTemplateOuterPos(currentIdx, fracStep);
  } else {
    // --- PHASE 4: STATIONARY AT REPLICATION STATION ---
    phi29Pos = { x: stationPos.x, y: stationPos.y };
  }

  // 5. TARGET STRAND DRAWING & PROGRESSIVE BASE-BY-BASE STRAND DISPLACEMENT
  if (p < 0.60) {
    let targetAlpha = 1.0;
    let targetDriftY = 0;
    if (p >= 0.42) {
      const u = (p - 0.42) / 0.14;
      targetAlpha = Math.max(0, 1.0 - u);
      targetDriftY = u * 90;
    }

    if (targetAlpha > 0.01) {
      ctx.save();
      ctx.globalAlpha = targetAlpha;

      const dispTargetPts = [];
      const isBaseDisplacedList = [];

      for (let i = 0; i < geoHs.tLen; i++) {
        let pairedPlpIdx = null;
        let stepsToReach = Infinity;

        if (state.leftTargetIdx !== -1 && i >= state.leftTargetIdx && i < state.leftTargetIdx + lLen) {
          const offset = i - state.leftTargetIdx;
          const k = (lLen - 1) - offset;
          pairedPlpIdx = k;
          stepsToReach = dockStartIdx - k;
        } else if (state.rightTargetIdx !== -1 && i >= state.rightTargetIdx && i < state.rightTargetIdx + rLen) {
          const offset = i - state.rightTargetIdx;
          const j = (rLen - 1) - offset;
          pairedPlpIdx = (total - rLen) + j;
          stepsToReach = (dockStartIdx + 1) + (rLen - 1 - j);
        }

        let isDisplaced = false;
        let peelDist = 0;

        if (p < 0.10) {
          isDisplaced = false;
        } else if (p >= 0.42) {
          isDisplaced = true;
          peelDist = 35 + targetDriftY;
        } else if (pairedPlpIdx !== null) {
          if (stepsTraveled >= stepsToReach) {
            isDisplaced = true;
            const stepsPast = stepsTraveled - stepsToReach;
            peelDist = Math.min(1.0, (stepsPast + 0.5) / 2.5) * 32;
          }
        } else {
          const distToLeft = (state.leftTargetIdx !== -1 && i < state.leftTargetIdx) ? (state.leftTargetIdx - i) : 999;
          const distToRight = (state.rightTargetIdx !== -1 && i >= state.rightTargetIdx + rLen) ? (i - (state.rightTargetIdx + rLen - 1)) : 999;
          if (distToLeft < 999) {
            const stepsToReachLeft = dockStartIdx - (lLen - 1);
            if (stepsTraveled >= stepsToReachLeft) {
              isDisplaced = true;
              peelDist = Math.min(1.0, (stepsTraveled - stepsToReachLeft) / 2.5) * 32;
            }
          } else if (distToRight < 999) {
            const stepsToReachRight = (dockStartIdx + 1) + (rLen - 1);
            if (stepsTraveled >= stepsToReachRight) {
              isDisplaced = true;
              peelDist = Math.min(1.0, (stepsTraveled - stepsToReachRight) / 2.5) * 32;
            }
          }
        }

        isBaseDisplacedList.push(isDisplaced);
        const wave = isDisplaced ? Math.sin(i * 0.4 + performance.now() * 0.005) * 3 : 0;
        dispTargetPts.push({
          x: geoHs.targetPts[i].x,
          y: geoHs.targetY + peelDist + targetDriftY + wave
        });
      }

      // Draw Target Strand with Sequence
      drawStrandWithSequence(dispTargetPts, state.target, "#f59e0b", getCanvasText('strands', 'targetStrandLabel', "Target ({targetType}) 5' → 3'", { targetType: state.targetType }), false);

      // Draw Watson-Crick Base-Pairing Rungs strictly for non-displaced bases
      if (isBasesVisible) {
        ctx.save();
        ctx.lineWidth = 1.8;
        ctx.setLineDash([2, 3]);

        // Left Arm Rungs
        for (let k = 0; k < lLen; k++) {
          const plpIdx = k;
          const tIdx = state.leftTargetIdx !== -1 ? (state.leftTargetIdx + (lLen - 1 - k)) : -1;
          if (tIdx >= 0 && tIdx < geoHs.tLen && !isBaseDisplacedList[tIdx]) {
            const armBase = (state.leftArmSeq && state.leftArmSeq[k]) || (state.plp && state.plp[plpIdx]) || 'N';
            const targetBase = state.target[tIdx] || 'N';
            const matched = isWatsonCrick(armBase, targetBase);
            ctx.strokeStyle = matched ? "rgba(52, 211, 153, 0.9)" : "rgba(244, 63, 94, 0.9)";
            ctx.beginPath();
            ctx.moveTo(plpCurrentPts[plpIdx].x, plpCurrentPts[plpIdx].y);
            ctx.lineTo(dispTargetPts[tIdx].x, dispTargetPts[tIdx].y);
            ctx.stroke();
          }
        }

        // Right Arm Rungs
        for (let j = 0; j < rLen; j++) {
          const plpIdx = (total - rLen) + j;
          const tIdx = state.rightTargetIdx !== -1 ? (state.rightTargetIdx + (rLen - 1 - j)) : -1;
          if (tIdx >= 0 && tIdx < geoHs.tLen && !isBaseDisplacedList[tIdx]) {
            const armBase = (state.rightArmSeq && state.rightArmSeq[j]) || (state.plp && state.plp[plpIdx]) || 'N';
            const targetBase = state.target[tIdx] || 'N';
            const matched = isWatsonCrick(armBase, targetBase);
            ctx.strokeStyle = matched ? "rgba(52, 211, 153, 0.9)" : "rgba(244, 63, 94, 0.9)";
            ctx.beginPath();
            ctx.moveTo(plpCurrentPts[plpIdx].x, plpCurrentPts[plpIdx].y);
            ctx.lineTo(dispTargetPts[tIdx].x, dispTargetPts[tIdx].y);
            ctx.stroke();
          }
        }
        ctx.restore();
      }

      // Status annotation
      if (p >= 0.10 && p < 0.42) {
        ctx.save();
        ctx.fillStyle = "#fbbf24";
        ctx.font = "bold 10px JetBrains Mono";
        ctx.textAlign = "center";
        ctx.fillText("⚡ PHI29 DISPLACING TARGET STRAND BASE-BY-BASE", geoHs.cx, geoHs.targetY + 50);
        ctx.restore();
      } else if (p >= 0.42 && targetAlpha > 0.1) {
        ctx.save();
        ctx.fillStyle = "#fbbf24";
        ctx.font = "bold 10px JetBrains Mono";
        ctx.textAlign = "center";
        ctx.fillText("✓ TARGET STRAND FULLY DISPLACED", geoHs.cx, geoHs.targetY + 70 + targetDriftY);
        ctx.restore();
      }

      ctx.restore();
    }
  }

  // 6. RCA PRIMER: ANNEALED WITH INTACT HYDROGEN BONDS (PHASES 1-3 ONLY)
  // In Phase 4 (p >= 0.56), the primer is peeled into the concatemer tail in Section 9!
  if (p < 0.56) {
    // A. Base-pairing rungs between primer and PLP template
    if (isBasesVisible && state.primer && state.plp) {
      ctx.save();
      ctx.lineWidth = 1.5;
      ctx.setLineDash([2, 3]);
      for (let k = 0; k < primerLen; k++) {
        const pB = state.primer[k];
        const tIdx = (dockStartIdx + (primerLen - 1 - k)) % total;
        const tB = state.plp[tIdx];
        const isMatch = isWatsonCrick(pB, tB);
        ctx.strokeStyle = isMatch ? "#10b981" : "#ef4444";
        ctx.beginPath();
        ctx.moveTo(primerPts[k].x, primerPts[k].y);
        ctx.lineTo(templatePts[k].x, templatePts[k].y);
        ctx.stroke();
      }
      ctx.restore();
    }

    // B. Primer strand drawing
    drawStrandWithSequence(primerPts, state.primer, "#f43f5e", "", false);

    // B2. Biotin badge at 5' terminus (primerPts[0])
    const bPt = primerPts[0];
    drawBiotinTag(ctx, bPt, bPt.x + 20, bPt.y - 18, 1.0, "Biotin");

    if (p >= 0.42 && p < 0.56) {
      ctx.save();
      ctx.fillStyle = "rgba(59, 130, 246, 0.25)";
      ctx.strokeStyle = "#3b82f6";
      ctx.lineWidth = 1;
      const bW = 340;
      const bH = 20;
      const bX = phi29Pos.x - bW * 0.5;
      const bY = phi29Pos.y + 36;
      ctx.beginPath();
      ctx.roundRect(bX, bY, bW, bH, 6);
      ctx.fill();
      ctx.stroke();

      ctx.fillStyle = "#ffffff";
      ctx.font = "bold 9px 'JetBrains Mono', monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("⚡ TARGET DISPLACED • TEMPLATE MORPHING TO REPLICATION STATION", phi29Pos.x, bY + bH * 0.5);
      ctx.restore();
    }
  }

  // 7. SYNTHESIZED DAUGHTER STRAND & DOUBLE-STRANDED CIRCLE (PHASE 2 & PHASE 3 & PHASE 4)
  let currentActiveDntp = null;
  let latestSynthPt = null;
  let totalDntpsAdded = 0;
  let closestIdx = 0;

  if (p >= 0.10 && p < 0.56) {
    // --- PHASE 2 & 3: SYNTHESIZING DAUGHTER STRAND EXTENDING FROM PRIMER 3'-OH ---
    const totalSteps = Math.floor(stepsTraveled);
    totalDntpsAdded = totalSteps;

    if (totalSteps > 0) {
      const synthPts = [];
      const synthBases = [];
      const tempPts = [];

      // Primer 3'-OH is at primerPts[primerLen - 1]
      const p3Pt = primerPts[primerLen - 1];

      for (let s = 1; s <= totalSteps; s++) {
        const tIdx = (dockStartIdx - s + total * 10) % total;
        const tPt = plpCurrentPts[tIdx];
        tempPts.push(tPt);

        let sx, sy;
        if (m < 0.5 && tPt.y >= geoHs.landingY - 14) {
          // Bottom arms: daughter strand is synthesized on outer/lower face (where target was displaced)
          sx = tPt.x;
          sy = tPt.y + baseDist;
        } else {
          const dx = tPt.x - currentCenter.x;
          const dy = tPt.y - currentCenter.y;
          const len = Math.hypot(dx, dy) || 1;
          sx = tPt.x + (dx / len) * baseDist;
          sy = tPt.y + (dy / len) * baseDist;
        }
        synthPts.push({ x: sx, y: sy });

        const tBase = (state.plp && state.plp[tIdx]) || 'A';
        const dBase = complementBase(tBase);
        synthBases.push(dBase);
      }

      // A. Synthesized strand backbone in vibrant cyan #38bdf8
      ctx.save();
      ctx.strokeStyle = isBasesVisible ? "rgba(56, 189, 248, 0.55)" : "#38bdf8";
      ctx.lineWidth = isBasesVisible ? 3.0 : 4.5;
      ctx.lineCap = "round";
      ctx.lineJoin = "round";
      ctx.shadowColor = "#38bdf8";
      ctx.shadowBlur = 10;
      ctx.beginPath();
      ctx.moveTo(p3Pt.x, p3Pt.y); // Extends directly from primer 3'-OH!
      for (let s = 0; s < synthPts.length; s++) {
        ctx.lineTo(synthPts[s].x, synthPts[s].y);
      }
      ctx.stroke();
      ctx.restore();

      // B. Watson-Crick hydrogen-bond rungs connecting template to newly incorporated dNTPs
      ctx.save();
      ctx.strokeStyle = "rgba(16, 185, 129, 0.85)";
      ctx.lineWidth = 1.6;
      ctx.setLineDash([2, 3]);
      for (let s = 0; s < synthPts.length; s++) {
        ctx.beginPath();
        ctx.moveTo(synthPts[s].x, synthPts[s].y);
        ctx.lineTo(tempPts[s].x, tempPts[s].y);
        ctx.stroke();
      }
      ctx.restore();

      // C. Render EVERY incorporated dNTP base along the daughter strand!
      for (let s = 0; s < synthPts.length; s++) {
        const baseLetter = synthBases[s];
        if (isBasesVisible) {
          drawNucleotideNode(ctx, synthPts[s].x, synthPts[s].y, baseLetter, currentRadius);
        } else {
          const bStyle = getBaseStyle(baseLetter);
          ctx.save();
          ctx.fillStyle = bStyle.bg;
          ctx.strokeStyle = bStyle.border;
          ctx.lineWidth = 1.2;
          ctx.beginPath();
          ctx.arc(synthPts[s].x, synthPts[s].y, 3.8, 0, Math.PI * 2);
          ctx.fill();
          ctx.stroke();
          ctx.restore();
        }
      }

      // D. Leading 3'-OH active incorporation site at Phi29
      if (synthPts.length > 0) {
        latestSynthPt = synthPts[synthPts.length - 1];
        const latestBase = synthBases[synthBases.length - 1];
        currentActiveDntp = "d" + latestBase + "TP";
        const bStyle = getBaseStyle(latestBase);

        const pulse = (Math.sin(performance.now() * 0.012) + 1) * 0.5;
        ctx.save();
        ctx.strokeStyle = bStyle.border;
        ctx.lineWidth = 2.5;
        ctx.shadowColor = bStyle.border;
        ctx.shadowBlur = 12 + pulse * 6;
        ctx.beginPath();
        ctx.arc(latestSynthPt.x, latestSynthPt.y, (isBasesVisible ? currentRadius : 4) + 3 + pulse * 2, 0, Math.PI * 2);
        ctx.stroke();
        ctx.restore();
      }
    }
  } else if (p >= 0.56) {
    // --- PHASE 4: DOUBLE-STRANDED CIRCULAR DNA WITH ACTIVE ROLLING REPLICATION ---
    // A. Active template base currently passing through phi29 active site at station
    closestIdx = 0;
    let minDiff = 999;
    for (let i = 0; i < total; i++) {
      let diff = Math.abs((geoCirc.plpPts[i].angle - stationAngle) % (Math.PI * 2));
      if (diff > Math.PI) diff = 2 * Math.PI - diff;
      if (diff < minDiff) {
        minDiff = diff;
        closestIdx = i;
      }
    }
    const activeTemplateBase = (state.plp && state.plp[closestIdx]) || 'A';
    const activeDntpBase = complementBase(activeTemplateBase);
    currentActiveDntp = "d" + activeDntpBase + "TP";
    const bStyle = getBaseStyle(activeDntpBase);

    // B. Calculate concentric outer ring points for the newly synthesized strand
    const outerPts = [];
    const compBases = [];
    for (let i = 0; i < total; i++) {
      const tPt = geoCirc.plpPts[i];
      const ang = tPt.angle;
      const ox = geoCirc.cx + (geoCirc.radius + baseDist) * Math.cos(ang);
      const oy = geoCirc.cy + (geoCirc.radius + baseDist) * Math.sin(ang);
      outerPts.push({ x: ox, y: oy, angle: ang });

      const tBase = (state.plp && state.plp[i]) || 'A';
      compBases.push(complementBase(tBase));
    }

    // C. Watson-Crick hydrogen-bond rungs (connecting template inner base to synthesized outer base)
    ctx.save();
    ctx.strokeStyle = "rgba(16, 185, 129, 0.85)";
    ctx.lineWidth = 1.6;
    ctx.setLineDash([2, 3]);
    for (let i = 0; i < total; i++) {
      ctx.beginPath();
      ctx.moveTo(geoCirc.plpPts[i].x, geoCirc.plpPts[i].y);
      ctx.lineTo(outerPts[i].x, outerPts[i].y);
      ctx.stroke();
    }
    ctx.restore();

    // D. Draw outer synthesized strand backbone (with a breakage at the active site)
    ctx.save();
    ctx.strokeStyle = isBasesVisible ? "rgba(56, 189, 248, 0.55)" : "#38bdf8";
    ctx.lineWidth = isBasesVisible ? 3.0 : 4.5;
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.shadowColor = "#38bdf8";
    ctx.shadowBlur = 10;
    ctx.beginPath();
    
    // Create a visual breakage (gap) between the newly synthesized 3' base and the displaced 5' base
    const displacedIdx = (closestIdx + 1) % total;
    const startArcAngle = outerPts[displacedIdx].angle;
    ctx.arc(geoCirc.cx, geoCirc.cy, geoCirc.radius + baseDist, startArcAngle, stationAngle + 2 * Math.PI, false);
    
    ctx.stroke();
    ctx.restore();

    // E. Draw nucleotide nodes on the outer synthesized ring
    for (let i = 0; i < total; i++) {
      const baseLetter = compBases[i];
      if (isBasesVisible) {
        drawNucleotideNode(ctx, outerPts[i].x, outerPts[i].y, baseLetter, currentRadius);
      } else {
        const bSt = getBaseStyle(baseLetter);
        ctx.save();
        ctx.fillStyle = bSt.bg;
        ctx.strokeStyle = bSt.border;
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.arc(outerPts[i].x, outerPts[i].y, 3.8, 0, Math.PI * 2);
        ctx.fill();
        ctx.stroke();
        ctx.restore();
      }
    }

    // F. Active incorporation site reference
    // (Rendered in foreground with Active Site Window in Section 8)
  }

  // 8. TRANSLUCENT GLOWING PHI29 ENZYME & FOREGROUND ACTIVE SITE WINDOW
  const pulse = (Math.sin(performance.now() * 0.012) + 1) * 0.5;

  // Compute orientation angle of DNA strand at phi29
  let normAngle = 0;
  if (p >= 0.56) {
    normAngle = Math.atan2(phi29Pos.y - geoCirc.cy, phi29Pos.x - geoCirc.cx);
  } else if (p >= 0.10) {
    normAngle = Math.atan2(phi29Pos.y - currentCenter.y, phi29Pos.x - currentCenter.x);
  }

  // --- A. TRANSLUCENT GLOWING ENZYME BODY ---
  // The crystal-clear glowing body allows the template strand, synthesized strand,
  // and base pairing underneath to be 100% visible!
  ctx.save();
  // Outer ambient halo
  ctx.fillStyle = "rgba(59, 130, 246, 0.20)";
  ctx.beginPath();
  ctx.arc(phi29Pos.x, phi29Pos.y, 30, 0, Math.PI * 2);
  ctx.fill();

  // Translucent glowing enzyme sphere
  const enzymeGrad = ctx.createRadialGradient(
    phi29Pos.x - 5, phi29Pos.y - 5, 2,
    phi29Pos.x, phi29Pos.y, 22
  );
  enzymeGrad.addColorStop(0, "rgba(96, 165, 250, 0.42)");
  enzymeGrad.addColorStop(0.65, "rgba(37, 99, 235, 0.32)");
  enzymeGrad.addColorStop(1.0, "rgba(29, 78, 216, 0.22)");

  ctx.fillStyle = enzymeGrad;
  ctx.strokeStyle = "#60a5fa";
  ctx.lineWidth = 2.2;
  ctx.shadowColor = "#3b82f6";
  ctx.shadowBlur = 18;
  ctx.beginPath();
  ctx.arc(phi29Pos.x, phi29Pos.y, 22, 0, Math.PI * 2);
  ctx.fill();
  ctx.stroke();

  // Subtle inner catalytic cleft ring
  ctx.strokeStyle = "rgba(147, 197, 253, 0.60)";
  ctx.lineWidth = 1.2;
  ctx.setLineDash([2, 3]);
  ctx.beginPath();
  ctx.arc(phi29Pos.x, phi29Pos.y, 11, 0, Math.PI * 2);
  ctx.stroke();

  // Energetic sparks pulsing from enzyme perimeter
  if (p >= 0.10) {
    const sparkPhase = performance.now() * 0.01;
    ctx.strokeStyle = "#93c5fd";
    ctx.lineWidth = 1.8;
    ctx.setLineDash([]);
    for (let a = 0; a < 4; a++) {
      const ang = sparkPhase + a * (Math.PI * 0.5);
      ctx.beginPath();
      ctx.moveTo(phi29Pos.x + Math.cos(ang) * 22, phi29Pos.y + Math.sin(ang) * 22);
      ctx.lineTo(phi29Pos.x + Math.cos(ang) * 29, phi29Pos.y + Math.sin(ang) * 29);
      ctx.stroke();
    }
  }

  ctx.restore();

  // --- B. FOREGROUND ACTIVE SITE WINDOW (OPTION 3) ---
  // Renders the template base & matching daughter base in crystal-clear foreground inside the groove
  if (p >= 0.56) {
    const activeTemplatePt = geoCirc.plpPts[closestIdx];
    const activeTemplateBase = (state.plp && state.plp[closestIdx]) || 'A';
    const activeDntpBase = complementBase(activeTemplateBase);
    const activeOuterX = geoCirc.cx + (geoCirc.radius + baseDist) * Math.cos(activeTemplatePt.angle);
    const activeOuterY = geoCirc.cy + (geoCirc.radius + baseDist) * Math.sin(activeTemplatePt.angle);

    const midX = (activeTemplatePt.x + activeOuterX) * 0.5;
    const midY = (activeTemplatePt.y + activeOuterY) * 0.5;

    // 1. Pulsing Catalytic Aura around the active base pair
    ctx.save();
    ctx.strokeStyle = `rgba(52, 211, 153, ${0.45 + pulse * 0.45})`;
    ctx.lineWidth = 1.8;
    ctx.setLineDash([3, 2]);
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 10 + pulse * 6;
    ctx.beginPath();
    ctx.arc(midX, midY, (isBasesVisible ? currentRadius : 4) + 11 + pulse * 2, 0, Math.PI * 2);
    ctx.stroke();
    ctx.restore();

    // 2. Active Watson-Crick hydrogen bonding rungs
    ctx.save();
    ctx.strokeStyle = "#10b981";
    ctx.lineWidth = 2.4;
    ctx.setLineDash([2, 2]);
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 8;
    ctx.beginPath();
    ctx.moveTo(activeTemplatePt.x, activeTemplatePt.y);
    ctx.lineTo(activeOuterX, activeOuterY);
    ctx.stroke();
    ctx.restore();

    // 3. Foreground template base node & synthesized dNTP node
    if (isBasesVisible) {
      drawNucleotideNode(ctx, activeTemplatePt.x, activeTemplatePt.y, activeTemplateBase, currentRadius, true);
      drawNucleotideNode(ctx, activeOuterX, activeOuterY, activeDntpBase, currentRadius, true);
    } else {
      drawNucleotideNode(ctx, activeTemplatePt.x, activeTemplatePt.y, activeTemplateBase, 4.2, true);
      drawNucleotideNode(ctx, activeOuterX, activeOuterY, activeDntpBase, 4.2, true);
    }

    // 4. Phosphodiester catalytic spark at the active 3'-OH incorporation site
    ctx.save();
    ctx.fillStyle = "#fef08a";
    ctx.shadowColor = "#eab308";
    ctx.shadowBlur = 10;
    ctx.beginPath();
    ctx.arc(activeOuterX, activeOuterY, 2.5 + pulse * 1.5, 0, Math.PI * 2);
    ctx.fill();
    ctx.restore();
  } else if (p >= 0.10 && p < 0.56 && latestSynthPt) {
    // Phase 2 & 3: Active base pair at Phi29 active cleft
    const latestTIdx = (dockStartIdx - totalDntpsAdded + total * 10) % total;
    const tPt = plpCurrentPts[latestTIdx];
    const tBase = (state.plp && state.plp[latestTIdx]) || 'A';
    const dBase = complementBase(tBase);

    // Active Watson-Crick bond
    ctx.save();
    ctx.strokeStyle = "#10b981";
    ctx.lineWidth = 2.2;
    ctx.setLineDash([2, 2]);
    ctx.shadowColor = "#10b981";
    ctx.shadowBlur = 8;
    ctx.beginPath();
    ctx.moveTo(tPt.x, tPt.y);
    ctx.lineTo(latestSynthPt.x, latestSynthPt.y);
    ctx.stroke();
    ctx.restore();

    if (isBasesVisible) {
      drawNucleotideNode(ctx, tPt.x, tPt.y, tBase, currentRadius, true);
      drawNucleotideNode(ctx, latestSynthPt.x, latestSynthPt.y, dBase, currentRadius, true);
    }
  }

  // Status label below phi29 in Phase 2
  if (p >= 0.10 && p < 0.42) {
    ctx.save();
    ctx.font = "bold 8px 'JetBrains Mono', monospace";
    ctx.textAlign = "center";
    ctx.fillStyle = "#67e8f9";
    ctx.shadowColor = "#06b6d4";
    ctx.shadowBlur = 6;
    ctx.fillText("⚡ Displacing Target Strand", phi29Pos.x, phi29Pos.y + 36);
    ctx.restore();
  }

  // 8.5 FLOATING dNTP PARTICLES IN REACTION SOLUTION & ACTIVE RECRUITMENT
  if (state.showDntps !== false) {
    updateAndDrawDntps(ctx, w, h, dt, phi29Pos, p >= 0.10 && state.isPlaying, currentActiveDntp);
  }

  // 9. CONTINUOUS CONCATEMER SPOOLING, DISPLACED PRIMER COLORING & PROBE HYBRIDIZATION
  if (p >= 0.56) {
    const rcaProgress = Math.min((p - 0.56) / 0.44, 1.0);
    const repeatBases = total; // Each monomer repeat is exactly 1 circular template length (74 nt)
    const repeatPitch = repeatBases * spacing;
    const maxConcatLength = userRepeats * repeatPitch;
    const concatemerLength = maxConcatLength * (p >= 1.0 ? 1.0 : rcaProgress);
    const numBases = Math.max(1, Math.floor(concatemerLength / spacing));
    const cPts = [];
    
    // Starts at the exact mathematical coordinate of the outer synthesized ring at stationAngle,
    // but pushed slightly away to show displacement breakage
    const gapOffset = spacing * 1.35; 
    const forkX = geoCirc.cx + (geoCirc.radius + baseDist) * Math.cos(stationAngle) + gapOffset;
    const forkY = geoCirc.cy + (geoCirc.radius + baseDist) * Math.sin(stationAngle);

    for (let i = 0; i <= numBases; i++) {
      const segX = forkX + i * spacing;
      // Damped wave envelope ensures concatemer is 100% rigidly anchored at fork (i=0)
      const waveEnvelope = Math.min(1.0, i / 6.0);
      const wave = Math.sin(i * 0.35 - state.rcaRotationAngle * 1.5) * 16 * waveEnvelope;
      cPts.push({ x: segX, y: forkY + wave });
    }

    // Concatemer ssDNA synthesized in the 3'->5' direction from circular PLP template
    // is the exact Watson-Crick complement of PLP:
    const compPlp = state.plp ? state.plp.split('').map(b => complementBase(b, false)).join('') : "A";

    // The displaced RCA primer is at the leading front (far right) of the growing strand
    const primerStartIdx = Math.max(0, numBases - primerLen);

    // RCA primer binding region along concatemer repeats:
    // From 5' tip (numBases), each repeat has length total; the primer/binding segment has length primerLen:
    const isPrimerRegion = (idx) => {
      const distFrom5 = (numBases - idx) % total;
      return distFrom5 >= 0 && distFrom5 < primerLen;
    };

    // Helper to get biochemically consistent, deterministic nucleotide along concatemer
    const getConcatBase = (idx) => {
      const distFrom5 = numBases - idx;
      if (distFrom5 < primerLen) {
        return (state.primer && state.primer[distFrom5]) || 'A';
      }
      const s = distFrom5 - primerLen + 1;
      const tIdx = (dockStartIdx - s + total * 100000) % total;
      return (state.plp && complementBase(state.plp[tIdx])) || 'A';
    };

    // --- DRAW STRAND BACKBONE (SEGMENTED: CONCATEMER CYAN #38bdf8 + RCA PRIMER REGIONS RED #f43f5e) ---
    ctx.save();
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    
    // 1. Draw connection from outer circle peeling point to concatemer start
    const inPrimerFirst = isPrimerRegion(0);
    ctx.strokeStyle = inPrimerFirst
      ? (isBasesVisible ? "rgba(244, 63, 94, 0.75)" : "#f43f5e")
      : (isBasesVisible ? "rgba(56, 189, 248, 0.45)" : "#38bdf8");
    ctx.lineWidth = inPrimerFirst ? (isBasesVisible ? 3.0 : 4.5) : (isBasesVisible ? 2.5 : 4);
    if (!isBasesVisible) {
      ctx.shadowColor = inPrimerFirst ? "#f43f5e" : "#38bdf8";
      ctx.shadowBlur = inPrimerFirst ? 12 : 8;
    } else {
      ctx.shadowBlur = inPrimerFirst ? 6 : 0;
    }
    
    // Connect exactly to the active site (yellow glowing base)
    const peelX = geoCirc.cx + (geoCirc.radius + baseDist) * Math.cos(stationAngle);
    const peelY = geoCirc.cy + (geoCirc.radius + baseDist) * Math.sin(stationAngle);
    
    ctx.beginPath();
    ctx.moveTo(peelX, peelY);
    ctx.lineTo(cPts[0].x, cPts[0].y);
    ctx.stroke();

    // 2. Draw rest of concatemer
    for (let i = 0; i < numBases; i++) {
      const inPrimer = isPrimerRegion(i + 1);
      ctx.strokeStyle = inPrimer
        ? (isBasesVisible ? "rgba(244, 63, 94, 0.75)" : "#f43f5e")
        : (isBasesVisible ? "rgba(56, 189, 248, 0.45)" : "#38bdf8");
      ctx.lineWidth = inPrimer ? (isBasesVisible ? 3.0 : 4.5) : (isBasesVisible ? 2.5 : 4);
      if (!isBasesVisible) {
        ctx.shadowColor = inPrimer ? "#f43f5e" : "#38bdf8";
        ctx.shadowBlur = inPrimer ? 12 : 8;
      } else {
        ctx.shadowBlur = inPrimer ? 6 : 0;
      }
      ctx.beginPath();
      ctx.moveTo(cPts[i].x, cPts[i].y);
      ctx.lineTo(cPts[i + 1].x, cPts[i + 1].y);
      ctx.stroke();
    }
    ctx.restore();

    // --- DRAW NUCLEOTIDE NODES ALONG STRAND (DYNAMIC CONTINUOUS SPOOLING) ---
    if (isBasesVisible) {
      const radius = renderSettings.radius;
      for (let i = 0; i <= numBases; i++) {
        const inPrimer = isPrimerRegion(i);
        const b = getConcatBase(i);
        const sub = (i === numBases ? "5'" : "");
        drawNucleotideNode(ctx, cPts[i].x, cPts[i].y, b, radius, inPrimer, sub);
      }
    } else {
      // Dots when zoomed out
      for (let i = 0; i <= numBases; i++) {
        const inPrimer = isPrimerRegion(i);
        if (i % 2 === 0 || inPrimer) {
          ctx.fillStyle = inPrimer ? "#f43f5e" : "#38bdf8";
          ctx.beginPath();
          ctx.arc(cPts[i].x, cPts[i].y, inPrimer ? 3 : 2.5, 0, Math.PI * 2);
          ctx.fill();
        }
      }
    }

    // Labels & Glowing Biotin reporter tag at concatemer tip
    if (numBases > 0) {
      const tipPt = cPts[numBases];
      drawBiotinTag(ctx, tipPt, tipPt.x + 24, tipPt.y - 14, 1.0, "Biotin");

      ctx.save();
      ctx.fillStyle = "#f43f5e";
      ctx.font = "bold 9px 'JetBrains Mono', monospace";
      ctx.textAlign = "right";
      ctx.fillText("5'-RCA Primer", tipPt.x - 8, tipPt.y - 12);
      ctx.restore();
    }
    if (cPts.length > 4) {
      ctx.save();
      ctx.fillStyle = "#38bdf8";
      ctx.font = "bold 10px JetBrains Mono";
      const midTailIdx = Math.min(Math.floor(numBases * 0.4), Math.max(1, primerStartIdx));
      const midTailPt = cPts[Math.max(1, midTailIdx)];
      ctx.fillText("Growing Concatemer", midTailPt.x, midTailPt.y - 18);
      ctx.restore();
    }

    // Repeat markers
    state.concatemerRepeats = Math.min(userRepeats, Math.floor(concatemerLength / repeatPitch) + 1);
    const displayRepeats = Math.min(userRepeats, Math.floor(concatemerLength / repeatPitch));
    const repCounter = document.getElementById('repeats-counter');
    if (repCounter) repCounter.innerText = `${displayRepeats} / ${userRepeats}`;

    for (let r = 1; r <= userRepeats; r++) {
      const markerX = forkX + r * repeatPitch;
      if (markerX <= forkX + concatemerLength) {
        ctx.save();
        ctx.strokeStyle = "#eab308";
        ctx.lineWidth = 1;
        ctx.setLineDash([2, 2]);
        ctx.beginPath();
        ctx.moveTo(markerX, forkY - 25);
        ctx.lineTo(markerX, forkY + 25);
        ctx.stroke();

        ctx.fillStyle = "#fde047";
        ctx.font = "bold 9px JetBrains Mono";
        ctx.textAlign = "center";
        ctx.fillText(getCanvasText('rollingCircle', 'repeatBadgePrefix', `Repeat ${r}`, { repeatNum: r }), markerX, forkY - 30);
        ctx.restore();
      }
    }

    // 10. DETECTION PROBE HYBRIDIZATION WITH EXACT BASE-PAIRING ALIGNMENT
    const shouldShowProbes = (state.protocol && (state.protocol.hasProbe || state.protocol.isProbeBound)) || state.currentStage === 6;
    if (shouldShowProbes && state.probe) {
      const probeLen = state.probe.length;

      // Build full concatemer sequence string for rigorous, exact subsequence matching
      let concatSeq = "";
      for (let i = 0; i <= numBases; i++) {
        concatSeq += getConcatBase(i);
      }

      const targetSeqSense = state.probe.split('').map(b => complementBase(b, false)).join('');
      const targetSeqAntisense = reverseComplement(state.probe, false);

      let matchedTargetSeq = targetSeqSense;
      let probeDisplaySeq = state.probe;
      if (concatSeq.indexOf(targetSeqSense) === -1 && concatSeq.indexOf(targetSeqAntisense) !== -1) {
        matchedTargetSeq = targetSeqAntisense;
        probeDisplaySeq = state.probe.split('').reverse().join('');
      }

      // Locate all exact complementary target cassettes along the concatemer
      const siteIndices = [];
      if (state.isProbeComplementary) {
        let searchPos = 0;
        while (searchPos <= concatSeq.length - probeLen) {
          const foundIdx = concatSeq.indexOf(matchedTargetSeq, searchPos);
          if (foundIdx === -1) break;
          siteIndices.push(foundIdx);
          searchPos = foundIdx + probeLen;
        }
      }

      let boundCount = 0;
      if (state.isProbeComplementary && siteIndices.length > 0) {
        for (let r = 0; r < siteIndices.length; r++) {
          const siteStartBaseIdx = siteIndices[r];
          if (siteStartBaseIdx + probeLen <= numBases) {
            boundCount++;
            const probePts = [];
            const probeConcatPts = [];
            for (let k = 0; k < probeLen; k++) {
              const cIdx = siteStartBaseIdx + k;
              const cPt = cPts[cIdx];
              probeConcatPts.push(cPt);
              probePts.push({
                x: cPt.x,
                y: cPt.y + (isBasesVisible ? renderSettings.radius * 2 + 14 : 16)
              });
            }

            // Complementary Probe hybridizes with verified 1:1 vertical Watson-Crick rungs
            drawStrandWithSequence(probePts, probeDisplaySeq, "#eab308", getCanvasText('detection', 'probeStrandLabel', "Detection Probe (5' → 3')"), false);

            // Base pairing rungs (perfect vertical alignment: probePts[k].x === probeConcatPts[k].x)
            if (isBasesVisible) {
              ctx.save();
              ctx.lineWidth = 1.6;
              ctx.setLineDash([2, 3]);
              ctx.strokeStyle = "#10b981"; // Emerald green Watson-Crick
              for (let k = 0; k < probeLen; k++) {
                ctx.beginPath();
                ctx.moveTo(probePts[k].x, probePts[k].y);
                ctx.lineTo(probeConcatPts[k].x, probeConcatPts[k].y);
                ctx.stroke();
              }
              ctx.restore();
            }

            // 5'-FAM Fluorophore tethered via flexible chemical linker (FAM does NOT base-pair to DNA)
            const famX = probePts[0].x - 18;
            const famY = probePts[0].y + (isBasesVisible ? renderSettings.radius + 16 : 18);

            ctx.save();
            ctx.strokeStyle = "#94a3b8";
            ctx.lineWidth = 1.5;
            ctx.setLineDash([2, 2]);
            ctx.beginPath();
            ctx.moveTo(probePts[0].x, probePts[0].y);
            ctx.lineTo(famX, famY);
            ctx.stroke();
            ctx.setLineDash([]);

            // Vibrant 5'-FAM fluorophore pulsing orb
            const pulse = Math.sin(performance.now() * 0.008 + r) * 3 + 7;
            ctx.fillStyle = "rgba(74, 222, 128, 0.4)";
            ctx.shadowColor = "#4ade80";
            ctx.shadowBlur = 20;
            ctx.beginPath();
            ctx.arc(famX, famY, pulse + 4, 0, Math.PI * 2);
            ctx.fill();

            const famGrad = ctx.createRadialGradient(famX - 2, famY - 2, 1, famX, famY, pulse);
            famGrad.addColorStop(0, "#bbf7d0");
            famGrad.addColorStop(0.5, "#4ade80");
            famGrad.addColorStop(1, "#16a34a");
            ctx.fillStyle = famGrad;
            ctx.strokeStyle = "#15803d";
            ctx.lineWidth = 1.2;
            ctx.beginPath();
            ctx.arc(famX, famY, pulse, 0, Math.PI * 2);
            ctx.fill();
            ctx.stroke();

            ctx.fillStyle = "#ffffff";
            ctx.font = "bold 8px 'JetBrains Mono', monospace";
            ctx.textAlign = "center";
            ctx.fillText(getCanvasText('detection', 'fluorophoreActiveBadge', "5'-FAM"), famX, famY + pulse + 12);
            ctx.restore();
          }
        }
      } else {
        // Non-complementary: floating unbound in solution below concatemer
        const numProbesToShow = Math.min(3, userRepeats);
        for (let r = 1; r <= numProbesToShow; r++) {
          const floatBaseX = forkX + r * (total * spacing * 0.7);
          const floatY = forkY + 70 + Math.sin(performance.now() * 0.002 + r * 1.5) * 10;
          const unboundPts = [];
          for (let k = 0; k < probeLen; k++) {
            unboundPts.push({ x: floatBaseX + k * spacing, y: floatY });
          }
          drawStrandWithSequence(unboundPts, state.probe, "#64748b", getCanvasText('detection', 'probeStrandLabel', "Unbound Probe"), false);

          const famX = unboundPts[0].x - 18;
          const famY = unboundPts[0].y + (isBasesVisible ? renderSettings.radius + 16 : 18);

          ctx.save();
          ctx.strokeStyle = "#64748b";
          ctx.lineWidth = 1.2;
          ctx.setLineDash([2, 2]);
          ctx.beginPath();
          ctx.moveTo(unboundPts[0].x, unboundPts[0].y);
          ctx.lineTo(famX, famY);
          ctx.stroke();
          ctx.setLineDash([]);

          ctx.fillStyle = "#475569";
          ctx.beginPath();
          ctx.arc(famX, famY, 6, 0, Math.PI * 2);
          ctx.fill();
          ctx.fillStyle = "#94a3b8";
          ctx.font = "bold 8px 'JetBrains Mono', monospace";
          ctx.textAlign = "center";
          ctx.fillText(getCanvasText('detection', 'fluorophoreQuenchedBadge', "FAM (Off)"), famX, famY + 14);
          ctx.restore();
        }
      }

      state.probesBound = boundCount;
      const fluoroValEl = document.getElementById('fluoro-val');
      if (state.isProbeComplementary) {
        state.fluoroSignal = Math.min(100, Math.floor(boundCount * (100 / userRepeats)));
        if (fluoroValEl) fluoroValEl.innerText = `${state.fluoroSignal} ${getCanvasText('detection', 'fluorometerIntensityUnit', 'a.u.')}`;
      } else {
        state.fluoroSignal = 0;
        if (fluoroValEl) fluoroValEl.innerText = `0 ${getCanvasText('detection', 'fluorometerIntensityUnit', 'a.u.')} (Quenched)`;
      }

      if (state.fluoroHistory.length > 50) state.fluoroHistory.shift();
      state.fluoroHistory.push(state.fluoroSignal);
    }
  }
}

// --- STAGE 6: DETECTION PROBE HYBRIDIZATION ---
function drawStage6(w, h, p, dt) {
  drawStage5(w, h, 1.0, dt);
}
