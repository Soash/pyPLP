/**
 * pyPLP Simulator - Nucleotide Node & Sequence Strand Rendering
 */

const BASE_STYLES = {
  'A': { bg: '#059669', border: '#34d399', text: '#ffffff' }, // Emerald
  'T': { bg: '#e11d48', border: '#fb7185', text: '#ffffff' }, // Rose/Red
  'U': { bg: '#e11d48', border: '#fb7185', text: '#ffffff' }, // Rose/Red
  'C': { bg: '#0284c7', border: '#38bdf8', text: '#ffffff' }, // Sky Blue
  'G': { bg: '#d97706', border: '#fcd34d', text: '#ffffff' }, // Amber
  'N': { bg: '#475569', border: '#94a3b8', text: '#ffffff' }
};

function getBaseStyle(char) {
  const c = (char || 'N').toUpperCase();
  return BASE_STYLES[c] || BASE_STYLES['N'];
}

function drawNucleotideNode(targetCtx, x, y, base, radius, isHighlighted = false, sublabel = '') {
  const c = targetCtx || ctx;
  const style = getBaseStyle(base);
  c.save();
  c.beginPath();
  c.arc(x, y, radius, 0, Math.PI * 2);
  c.fillStyle = style.bg;
  c.fill();
  c.lineWidth = isHighlighted ? 2.5 : 1.2;
  c.strokeStyle = isHighlighted ? '#ffffff' : style.border;
  c.stroke();

  // Base letter
  c.fillStyle = style.text;
  const fontSize = Math.max(8, Math.floor(radius * 1.15));
  c.font = `bold ${fontSize}px 'JetBrains Mono', monospace`;
  c.textAlign = 'center';
  c.textBaseline = 'middle';
  c.fillText(base.toUpperCase(), x, y);

  // Optional sublabel (e.g. 5', 3'-OH)
  if (sublabel) {
    c.font = `bold 8px 'JetBrains Mono', monospace`;
    c.fillStyle = isHighlighted ? '#38bdf8' : '#94a3b8';
    c.fillText(sublabel, x, y - radius - 6);
  }
  c.restore();
}

// Strand renderer supporting zoom & settings: smooth strand when zoomed out, individual bases when zoomed in
function drawStrandWithSequence(points, sequence, defaultColor, label = "", labelAtEnd = false) {
  if (points.length < 2) return;
  const c = ctx;
  c.save();

  const isBasesVisible = shouldDrawBases();

  // Draw backbone line
  c.strokeStyle = isBasesVisible ? "rgba(148, 163, 184, 0.45)" : defaultColor;
  c.lineWidth = isBasesVisible ? 2.5 : 4;
  c.lineCap = "round";
  c.lineJoin = "round";
  if (!isBasesVisible) {
    c.shadowColor = defaultColor;
    c.shadowBlur = 10;
  }

  c.beginPath();
  c.moveTo(points[0].x, points[0].y);
  for (let i = 1; i < points.length; i++) {
    c.lineTo(points[i].x, points[i].y);
  }
  c.stroke();
  c.shadowBlur = 0;

  // When zoomed in or always on, draw individual nucleotide bases
  if (isBasesVisible && sequence && sequence.length > 0) {
    const radius = renderSettings.radius;
    for (let i = 0; i < sequence.length; i++) {
      const t = i / Math.max(sequence.length - 1, 1);
      const ptIdx = Math.min(Math.floor(t * (points.length - 1)), points.length - 2);
      const localT = (t * (points.length - 1)) - ptIdx;
      const p1 = points[ptIdx];
      const p2 = points[ptIdx + 1];
      const x = p1.x + (p2.x - p1.x) * localT;
      const y = p1.y + (p2.y - p1.y) * localT;

      const base = sequence[i];
      let sub = '';
      if (i === 0) sub = getCanvasText('strands', 'nucleotide5PrimeTag', "5'");
      else if (i === sequence.length - 1) sub = getCanvasText('strands', 'nucleotide3PrimeTag', "3'");

      drawNucleotideNode(c, x, y, base, radius, false, sub);
    }
  } else if (!isBasesVisible) {
    // Just small ticks along strand
    c.fillStyle = defaultColor;
    for (let i = 0; i < points.length; i += 3) {
      c.beginPath();
      c.arc(points[i].x, points[i].y, 2.5, 0, Math.PI * 2);
      c.fill();
    }
  }

  // Draw label
  if (label) {
    c.fillStyle = defaultColor;
    c.font = "bold 11px JetBrains Mono";
    const lPt = labelAtEnd ? points[points.length - 1] : points[0];
    c.fillText(label, lPt.x - 20, lPt.y - 14);
  }
  c.restore();
}

function drawBackgroundGrid(w, h) {
  const c = ctx;
  c.strokeStyle = "rgba(51, 65, 85, 0.2)";
  c.lineWidth = 1;
  const step = 40;
  const startX = -1000;
  const endX = w + 1000;
  const startY = -1000;
  const endY = h + 1000;

  for (let x = startX; x < endX; x += step) {
    c.beginPath();
    c.moveTo(x, startY);
    c.lineTo(x, endY);
    c.stroke();
  }
  for (let y = startY; y < endY; y += step) {
    c.beginPath();
    c.moveTo(startX, y);
    c.lineTo(endX, y);
    c.stroke();
  }
}
