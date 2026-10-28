/**
 * pyPLP Simulator - Biophysical Analysis Engine
 */

function complementBase(b, isRNA = false) {
  b = (b || 'N').toUpperCase();
  if (b === 'A') return isRNA ? 'U' : 'T';
  if (b === 'T' || b === 'U') return 'A';
  if (b === 'G') return 'C';
  if (b === 'C') return 'G';
  return 'N';
}

function isWatsonCrick(b1, b2) {
  b1 = (b1 || '').toUpperCase();
  b2 = (b2 || '').toUpperCase();
  if ((b1 === 'A' && (b2 === 'T' || b2 === 'U')) || ((b1 === 'T' || b1 === 'U') && b2 === 'A')) return true;
  if ((b1 === 'G' && b2 === 'C') || (b1 === 'C' && b2 === 'G')) return true;
  return false;
}

function reverseComplement(seq, isRNA = false) {
  return seq.split('').reverse().map(b => complementBase(b, isRNA)).join('');
}

function estimateTm(seq) {
  if (!seq || seq.length === 0) return 0;
  let gc = 0;
  for (let char of seq.toUpperCase()) {
    if (char === 'G' || char === 'C') gc++;
  }
  const gcPct = (gc / seq.length) * 100;
  const tm = 81.5 + 16.6 * Math.log10(0.05 + 0.12 * 1.5) + 0.41 * gcPct - (675 / seq.length);
  return Math.round(tm * 10) / 10;
}

function findBestSubsequenceMatch(query, target) {
  if (!query || !target || query.length === 0 || target.length === 0) {
    return { found: false, matchCount: 0, total: 0, pct: 0, startIdx: -1 };
  }
  const q = query.toUpperCase();
  const t = target.toUpperCase();
  const qLen = q.length;
  if (t.length < qLen) {
    return { found: false, matchCount: 0, total: qLen, pct: 0, startIdx: -1 };
  }
  let bestMatches = 0;
  let bestIdx = -1;
  for (let i = 0; i <= t.length - qLen; i++) {
    let matches = 0;
    for (let j = 0; j < qLen; j++) {
      if (q[j] === t[i + j]) matches++;
    }
    if (matches > bestMatches) {
      bestMatches = matches;
      bestIdx = i;
    }
  }
  const pct = Math.round((bestMatches / qLen) * 100);
  return {
    found: (bestMatches / qLen) >= 0.70, // 70% threshold required to initiate annealing
    matchCount: bestMatches,
    total: qLen,
    pct: pct,
    startIdx: bestIdx
  };
}

function analyzeMolecules() {
  const plp = state.plp.toUpperCase();
  const target = state.target.toUpperCase().replace(/U/g, state.targetType === 'RNA' ? 'U' : 'T');
  const primer = state.primer.toUpperCase();
  const probe = state.probe.toUpperCase();
  
  state.leftArmLen = 15;
  state.rightArmLen = 15;
  if (plp.length >= 30) {
    state.leftArmSeq = plp.slice(0, 15);
    state.rightArmSeq = plp.slice(plp.length - 15);
    state.backboneSeq = plp.slice(15, plp.length - 15);
  } else {
    state.leftArmSeq = plp;
    state.rightArmSeq = "";
    state.backboneSeq = "";
  }

  const targetRevComp = reverseComplement(target, false);
  const leftArmRC = reverseComplement(state.leftArmSeq, false);
  const rightArmRC = reverseComplement(state.rightArmSeq, false);
  
  let leftMatch = target.indexOf(leftArmRC);
  let rightMatch = target.indexOf(rightArmRC);

  // If exact match not found, check best match (e.g. for mismatch preset)
  if (leftMatch === -1 && state.leftArmSeq.length >= 10) {
    const bestL = findBestSubsequenceMatch(leftArmRC, target);
    if (bestL.found) leftMatch = bestL.startIdx;
  }
  if (rightMatch === -1 && state.rightArmSeq.length >= 10) {
    const bestR = findBestSubsequenceMatch(rightArmRC, target);
    if (bestR.found) rightMatch = bestR.startIdx;
  }

  state.leftTargetIdx = leftMatch;
  state.rightTargetIdx = rightMatch;

  let hasNickMismatch = false;
  if (leftMatch !== -1 && rightMatch !== -1) {
    const expectedRightStart = leftMatch + state.leftArmSeq.length;
    state.nickGap = rightMatch - expectedRightStart;

    // Check Watson-Crick pairing at the critical nick junction:
    // Left Arm 5'-PO4 (base 0) pairs with target at (leftMatch + leftArmSeq.length - 1)
    const leftJunctionTarget = (leftMatch + state.leftArmSeq.length - 1 < target.length) ? target[leftMatch + state.leftArmSeq.length - 1] : 'N';
    const leftJunctionMatch = isWatsonCrick(state.leftArmSeq[0], leftJunctionTarget);

    // Right Arm 3'-OH (last base) pairs with target at rightMatch
    const rightJunctionTarget = (rightMatch < target.length) ? target[rightMatch] : 'N';
    const rightJunctionMatch = isWatsonCrick(state.rightArmSeq[state.rightArmSeq.length - 1], rightJunctionTarget);

    if (!leftJunctionMatch || !rightJunctionMatch) {
      hasNickMismatch = true;
    }

    state.isLigatable = (state.nickGap === 0 && !hasNickMismatch);
    state.hasNickMismatch = hasNickMismatch;
  } else {
    state.nickGap = 999;
    state.isLigatable = false;
    state.hasNickMismatch = false;
  }

  state.leftTm = estimateTm(state.leftArmSeq);
  state.rightTm = estimateTm(state.rightArmSeq);
  state.deltaTm = Math.abs(Math.round((state.leftTm - state.rightTm) * 10) / 10);

  // --- RCA PRIMER COMPLEMENTARITY ANALYSIS ---
  // Circularized padlock template: primer (5'->3') must anneal ANTIPARALLEL to circular template (5'->3').
  // Therefore, the target binding cassette on the circular template must be reverseComplement(primer).
  const circPlp = plp + plp;
  const primerRC = reverseComplement(primer, false);
  const bestPrimerMatch = findBestSubsequenceMatch(primerRC, circPlp);

  state.primerTm = estimateTm(primer);
  if (bestPrimerMatch.found && bestPrimerMatch.pct >= 70) {
    state.isPrimerComplementary = true;
    state.primerMatchPct = bestPrimerMatch.pct;
    state.primerDockIdx = bestPrimerMatch.startIdx % Math.max(plp.length, 1);
  } else {
    state.isPrimerComplementary = false;
    state.primerMatchPct = bestPrimerMatch.pct;
    state.primerDockIdx = -1;
  }

  // --- DETECTION PROBE COMPLEMENTARITY ANALYSIS ---
  // Concatemer is tandem repeats of reverseComplement(plp)
  // A probe hybridizes if its RC matches the concatemer (meaning probe matches plp)
  const probeMatchSense = findBestSubsequenceMatch(probe, circPlp);
  const probeMatchAntisense = findBestSubsequenceMatch(reverseComplement(probe, false), circPlp);
  const bestProbeMatch = (probeMatchSense.pct >= probeMatchAntisense.pct) ? probeMatchSense : probeMatchAntisense;

  state.probeTm = estimateTm(probe);
  if (bestProbeMatch.found) {
    state.isProbeComplementary = true;
    state.probeMatchPct = bestProbeMatch.pct;
    state.probeDockIdx = bestProbeMatch.startIdx % Math.max(plp.length, 1);
  } else {
    state.isProbeComplementary = false;
    state.probeMatchPct = bestProbeMatch.pct;
    state.probeDockIdx = -1;
  }

  if (typeof updateHUD === 'function') {
    updateHUD();
  }
}
