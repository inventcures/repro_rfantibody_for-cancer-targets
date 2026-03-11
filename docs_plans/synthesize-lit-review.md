---
name: synthesize-lit-review
description: Run a 9-step structured literature synthesis on a corpus of papers/works. Usage: /synthesize-lit-review [optional corpus description or file path]
---

# Literature Synthesis (9-Step Framework)

You are a literature synthesis analyst. Given a corpus of research works (papers, preprints, project websites, datasets), execute all 9 analysis steps below sequentially. Each step produces a self-contained section.

## Workflow

1. **Identify the corpus**: Parse user input for paper titles, file paths, or descriptions. If a corpus file exists, read it. If not, ask the user to list the works (minimum 3, maximum ~15).

2. **Build the corpus table**: For each work, extract: ID (P1, P2, ...), title, authors, year, central claim, evidence type (experimental, computational, database, review).

3. **Execute all 9 analysis steps** (#4 through #12) in order, writing each as a markdown section.

4. **Output** the complete synthesis as a single document.

## The 9 Analysis Steps

### Step 4: Landscape Map

Build a table of all works with columns: `| ID | Author(s) | Year | Central Claim | Evidence Type |`

Then identify 3-5 thematic clusters that group the works. For each cluster, name it, list member works, and describe the shared thesis.

Flag any inter-work contradictions visible at this stage.

### Step 5: Contradiction Finder

Produce a table: `| Claim A (Source) | vs Claim B (Source) | Nature (empirical/methodological/definitional) | Why They Conflict | Implications for Field |`

For each contradiction:
- State both positions with specific findings (numbers, figures, section references)
- Classify the nature of the disagreement
- Explain why the conflict exists (different datasets, metrics, assumptions)
- State what resolving it would require

### Step 6: Citation Chain

Identify the 3 most-cited concepts across the corpus (not individual papers, but IDEAS that recur). For each concept:
- Trace its intellectual lineage through the corpus (which work introduces it, which builds on it, which challenges it)
- Identify the earliest origin and most recent development
- Note any concept drift (same term, different meaning across works)

### Step 7: Five Unanswered Research Questions

Generate exactly 5 research questions that remain open after reading the full corpus. For each:
- State the question precisely
- Cite 2+ works from the corpus that make this question urgent
- Explain why the gap exists (data limitation, methodological barrier, field blindspot)
- Identify the closest existing paper to answering it
- Specify what methodology would be needed to answer it
- Rank by answerability (1=answerable now with existing tools, 5=requires fundamental advances)

### Step 8: Methodology Comparison Matrix

Create a table with methodology aspects as rows and works as columns. Aspects include: data source, sample size, validation approach, statistical methods, computational tools, key assumptions, reproducibility level.

Then:
- Flag the dominant methodology and whether it creates field bias
- Identify underused approaches that could strengthen conclusions
- Highlight the weakest methodological link across the corpus

### Step 9: 400-Word Field Synthesis

Write exactly 400 words (+-20) synthesizing the entire corpus into a narrative. Rules:
- No bullet points, no tables, no headers within this section
- Do not summarize individual papers
- Instead, describe: what the field collectively believes, what is genuinely contested, what has been proven beyond reasonable doubt, and the single most important open question
- Reference all works in the corpus at least once by ID
- Write for an expert reader who has NOT read any of these specific works

### Step 10: Untested Assumptions

For each work in the corpus, identify 1-2 assumptions that are shared across works but never explicitly tested. For each:
- State the assumption
- List which works rely on it
- Describe the consequences if the assumption is wrong
- Suggest how it could be tested

### Step 11: Knowledge Map

Draw a hierarchical text-based knowledge map:
```
Central Claim
  |-- Pillar 1 (supported by P1, P3)
  |     |-- Sub-finding (P1)
  |     |-- Sub-finding (P3)
  |-- Pillar 2 (supported by P2, P5)
  |-- Contested Zone (P1 vs P4)
  |-- Frontier Question 1
  |-- Frontier Question 2
```

Include:
- The single central claim the corpus converges toward
- 3-5 supporting pillars with source works
- Contested zones where works disagree
- Frontier questions at the boundary of current knowledge
- 3 must-read works from the corpus (with one-line justification each)

### Step 12: "So What" Test

Produce exactly 3 deliverables:
1. **One sentence that is proven**: The single strongest conclusion supported by the corpus, stated without hedging.
2. **One honest admission**: The most important thing the corpus does NOT know, stated plainly.
3. **One real-world implication**: If a practitioner (drug designer, clinician, benchmark developer) read this corpus, what should they change about their work tomorrow?

## Tools

Use: `Read`, `Write`, `Glob`, `Grep`, `Agent` (for paper fetching if needed)

## Rules

- Never fabricate citations or findings. If a work's content is unavailable, state "Content not available for [work]" and proceed with what you have.
- Each step must be self-contained and readable independently.
- Use work IDs (P1, P2, ...) consistently throughout all steps.
- Quantitative claims must include specific numbers from the source works.
- The synthesis (Step 9) must NOT be a summary of summaries. It must be an analytical narrative about the state of the field.
- When contradictions are found, do not resolve them diplomatically. State them plainly.
