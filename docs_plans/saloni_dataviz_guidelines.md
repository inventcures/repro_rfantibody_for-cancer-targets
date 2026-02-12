# Saloni's Guide to Data Visualization — Reference for All Report Figures
> Source: https://www.scientificdiscovery.dev/p/salonis-guide-to-data-visualization
> Applied to: ALL figures in the cancer targets tech report

---

## Mandatory Rules (violating any of these = redo the figure)

### Chart Type Selection
- Choose chart type that directly answers the specific question being asked
- Different metrics (absolute counts, rates, percentages) reveal different things — pick the right one
- Avoid chart types that obscure patterns (no pie charts, no 3D bar charts)

### Text & Labels
- **ALL text must be horizontal** — no rotated axis labels, no vertical text
- Use plain language for broad audiences; keep technical terms only when they distinguish between similar concepts
- Combine precision with accessibility: "average (median) value" style
- **Direct label** data elements instead of using separate legends (put category names next to the data)
- Use legends ONLY when categories repeat multiple times across the visualization
- Order categories logically (inherent sequence) or alphabetically — never randomly

### Color
- Match colors to familiar associations (red for problems/failures, green for success, blue for positive)
- **Must be colorblind-accessible** — test with Coblis or use established CB-safe palettes
- Use direct labeling to supplement color for CB viewers
- Remove colors when grayscale conveys the same information
- No colored text backgrounds that distract from data

### Axes & Scales
- Never display metrics without units
- Y-axis: include buffer space above max and below min data points
- Anchor near-zero data to zero; preserve spacing for data far from zero
- Never use extreme zoom that exaggerates trivial changes
- Never use full-scale (0-100%) when it minimizes meaningful variation
- **No double-axis charts** with mismatched scales (manufactures false correlations)

### Standalone Requirement
- Every chart must be self-contained: title, subtitle, units, data source ON the chart
- Do NOT relegate critical info to captions or body text only
- Charts will be shared without context — they must make sense alone

### Complexity Handling
- "Clearer" not "simpler" — complexity is fine if it aids understanding
- Use **small multiples (faceting)** to split crowded charts into panels
- Annotate and guide viewers through unfamiliar chart types
- Highlight key takeaways adjacent to the chart

### Things to NEVER Do
- No 3D rendering (especially bar charts) — complicates value reading
- No arrows on line charts implying future trajectory
- No confidence intervals when prediction intervals or raw data would be clearer
- No decorative elements competing with data
- No tilted, hyphenated, or rotated text anywhere

### Reproducibility
- Include data source attribution on every figure
- Document methodological choices
- Figures must be regenerable from the underlying data

---

## Application to Cancer Targets Tech Report

### Figure Types We'll Use

| Data | Chart Type | Saloni Principle |
|------|-----------|-----------------|
| Score distributions (pAE, RMSD, ddG) | Histograms with direct labels, small multiples per campaign | Faceting for clarity |
| Design funnel (total → filtered → ranked) | Horizontal bar or waterfall chart | Meaningful chart type for attrition |
| Cross-campaign pass rates | Horizontal bar chart, sorted by pass rate | Logical ordering |
| Per-stage timing breakdown | Stacked horizontal bar | Shows composition + total |
| Top candidates per campaign | Clean table (not chart) | Tables for exact values |
| Composite score vs individual metrics | Scatter plot with direct labels | Shows correlations |
| Target properties vs success rate | Scatter with annotations | Guides interpretation |

### Color Palette (CB-safe, consistent across all figures)

```python
# Colorblind-safe palette (Paul Tol's bright scheme)
COLORS = {
    'pass': '#228833',      # green — success
    'fail': '#CC3311',      # red — failure
    'rfdiffusion': '#0077BB',  # blue — stage 1
    'proteinmpnn': '#EE7733',  # orange — stage 2
    'rf2': '#AA3377',          # purple — stage 3
    'neutral': '#BBBBBB',     # gray — background/reference
}
```

### Matplotlib Style Settings

```python
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.family': 'Helvetica',
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
})
```

### 6 Guiding Questions (check every figure against these)
1. Is my chart type meaningful for this specific question?
2. Can I make it clearer (not simpler — clearer)?
3. If complicated, did I guide viewers through it?
4. Does it work standalone (title, subtitle, units, source)?
5. Is the presentation choice justifiable?
6. Is it reproducible from the data?
