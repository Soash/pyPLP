/**
 * pyPLP Simulator - Geometric Computation Engines (Unified Topology)
 */

function computeHorseshoeGeometry(w, h, p = 1.0) {
  const cx = w * 0.35; // Positioned to leave ample room for concatemer growth to the right
  const targetY = h * 0.67;
  const spacing = renderSettings.spacing;
  const isBasesVisible = shouldDrawBases();
  const radius = renderSettings.radius;

  const rLen = (state.rightArmSeq && state.rightArmSeq.length) || 15;
  const lLen = (state.leftArmSeq && state.leftArmSeq.length) || 15;
  const bbLen = (state.backboneSeq && state.backboneSeq.length) || 40;
  const tLen = (state.target && state.target.length) || 50;
  const totalPlpLen = lLen + bbLen + rLen;

  const nickX = cx;
  const gapOffset = state.isLigatable ? 0 : Math.max(state.nickGap === 999 ? 1 : state.nickGap, 0) * spacing;
  
  let targetStartX;
  if (state.leftTargetIdx !== -1) {
    const targetNickIndex = state.leftTargetIdx + lLen - 1;
    targetStartX = (nickX - spacing * 0.5) - (targetNickIndex * spacing);
  } else if (state.rightTargetIdx !== -1) {
    targetStartX = (nickX + spacing * 0.5) - (state.rightTargetIdx * spacing);
  } else {
    targetStartX = cx - ((tLen - 1) * spacing) * 0.5;
  }

  const targetPts = [];
  for (let i = 0; i < tLen; i++) {
    targetPts.push({ x: targetStartX + i * spacing, y: targetY });
  }

  const landingY = targetY - (isBasesVisible ? (radius * 2 + 18) : 28);
  const initialFloatY = h * 0.28;

  // Fully formed final horseshoe reference points at landing position
  const lArmFinal = [];
  for (let k = 0; k < lLen; k++) {
    const lx = (nickX - spacing * 0.5) - k * spacing;
    lArmFinal.push({ x: lx, y: landingY });
  }

  const rArmFinal = [];
  for (let j = 0; j < rLen; j++) {
    const rx = (nickX + spacing * 0.5 + gapOffset) + (rLen - 1 - j) * spacing;
    rArmFinal.push({ x: rx, y: landingY });
  }

  const bbStartFinal = lArmFinal[lLen - 1];
  const bbEndFinal = rArmFinal[0];
  const loopSpan = bbEndFinal.x - bbStartFinal.x;
  const loopHeightFinal = Math.min(230, Math.max(130, loopSpan * 0.38));

  const cp1Final = { x: bbStartFinal.x - 30, y: bbStartFinal.y - loopHeightFinal };
  const cp2Final = { x: bbEndFinal.x + 30, y: bbEndFinal.y - loopHeightFinal };

  const bbPtsFinal = [];
  for (let i = 0; i < bbLen; i++) {
    const t = bbLen <= 1 ? 0 : i / (bbLen - 1);
    const mt = 1 - t;
    const bx = mt*mt*mt * bbStartFinal.x + 3*mt*mt*t * cp1Final.x + 3*mt*t*t * cp2Final.x + t*t*t * bbEndFinal.x;
    const by = mt*mt*mt * bbStartFinal.y + 3*mt*mt*t * cp1Final.y + 3*mt*t*t * cp2Final.y + t*t*t * bbEndFinal.y;
    bbPtsFinal.push({ x: bx, y: by });
  }

  // Straight linear PLP reference points (matching Stage 1 linear strand layout)
  const straightStartX = cx - ((totalPlpLen - 1) * spacing) * 0.5;

  // Folding & Landing Timeline Mechanics:
  // Phase 1 (p in 0.0 -> 0.58): Straight PLP bends & loops into horseshoe shape floating above target
  // Phase 2 (p in 0.58 -> 1.0): Formed horseshoe descends onto target strand and anneals
  const pFold = Math.min(1.0, p / 0.58);
  const u = 0.5 - 0.5 * Math.cos(Math.PI * pFold); // S-curve folding progress (0 to 1)

  const pDesc = Math.max(0.0, Math.min(1.0, (p - 0.58) / 0.42));
  const wDesc = 0.5 - 0.5 * Math.cos(Math.PI * pDesc); // Descent progress (0 to 1)

  const currentArmY = initialFloatY + wDesc * (landingY - initialFloatY);
  const yShift = currentArmY - landingY;

  // Left Arm dynamic interpolation from straight line to folded arm
  const lArmPts = [];
  for (let k = 0; k < lLen; k++) {
    const strX = straightStartX + k * spacing;
    const strY = initialFloatY;
    const hsX = lArmFinal[k].x;
    const hsY = lArmFinal[k].y + yShift;
    lArmPts.push({
      x: (1 - u) * strX + u * hsX,
      y: (1 - u) * strY + u * hsY
    });
  }

  // Right Arm dynamic interpolation from straight line to folded arm
  const rArmPts = [];
  for (let j = 0; j < rLen; j++) {
    const strX = straightStartX + (lLen + bbLen + j) * spacing;
    const strY = initialFloatY;
    const hsX = rArmFinal[j].x;
    const hsY = rArmFinal[j].y + yShift;
    rArmPts.push({
      x: (1 - u) * strX + u * hsX,
      y: (1 - u) * strY + u * hsY
    });
  }

  // Backbone loop dynamic interpolation: arches upward from straight center into curved loop
  const bbPts = [];
  for (let i = 0; i < bbLen; i++) {
    const strX = straightStartX + (lLen + i) * spacing;
    const strY = initialFloatY;
    const hsX = bbPtsFinal[i].x;
    const hsY = bbPtsFinal[i].y + yShift;
    bbPts.push({
      x: (1 - u) * strX + u * hsX,
      y: (1 - u) * strY + u * hsY
    });
  }

  const bbStart = lArmPts[lLen - 1];
  const bbEnd = rArmPts[0];
  const loopHeight = loopHeightFinal * u;
  const cp1 = {
    x: bbStart.x - 30 * u,
    y: (1 - u) * initialFloatY + u * (cp1Final.y + yShift)
  };
  const cp2 = {
    x: bbEnd.x + 30 * u,
    y: (1 - u) * initialFloatY + u * (cp2Final.y + yShift)
  };

  // Full PLP points in order from 5' to 3'
  const plpPts = [];
  for (let i = 0; i < lLen; i++) plpPts.push(lArmPts[i]);
  for (let i = 0; i < bbLen; i++) plpPts.push(bbPts[i]);
  for (let i = 0; i < rLen; i++) plpPts.push(rArmPts[i]);

  return {
    cx, targetY, landingY, currentArmY, nickX, gapOffset, spacing, radius, isBasesVisible,
    lLen, bbLen, rLen, totalPlpLen, tLen,
    targetPts, lArmPts, rArmPts, bbPts, plpPts,
    bbStart, bbEnd, cp1, cp2,
    loopHeight
  };
}

function computeCircularGeometry(w, h, rotAngle = 0) {
  const cx = w * 0.35; // Shared X center with horseshoe
  const targetY = h * 0.67;
  const spacing = renderSettings.spacing;
  const total = (state.plp && state.plp.length) || 70;
  const defaultRadius = Math.max(90, Math.min(140, (total * spacing * 0.45) / (2 * Math.PI)));
  const radius = (state.stage5CircleRadius != null && !isNaN(state.stage5CircleRadius))
    ? Math.max(50, Math.min(260, Number(state.stage5CircleRadius)))
    : defaultRadius;
  const bRadius = Math.min(renderSettings.radius, Math.max(5, (2 * Math.PI * radius) / (total * 2.3)));
  const cy = (targetY - 35) - radius; // Bottom lands near landingY

  const plpPts = [];
  for (let i = 0; i < total; i++) {
    // Base 0 (5'-PO4) is at bottom (+PI/2) meeting Right Arm 3'-OH
    const angle = (i / total) * Math.PI * 2 + Math.PI * 0.5 + rotAngle;
    const x = cx + radius * Math.cos(angle);
    const y = cy + radius * Math.sin(angle);
    plpPts.push({ x, y, angle });
  }

  return { cx, cy, radius, bRadius, total, plpPts, exitX: cx, exitY: cy - radius };
}
