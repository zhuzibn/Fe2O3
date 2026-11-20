# Project Guidelines & Documentation

All contributors (including AI Agents) must follow these protocols. Record meaningful changes concisely under the relevant sections.

---

## Section 1: Source Code Changes (Changelog)
Update `CHANGELOG.md` for meaningful behavior changes or significant fixes,
following the global changelog policy. Group related edits into one entry.

Include meaningful features, breaking changes, and setup or migration impacts.
Skip typo fixes, formatting, behavior-preserving internal refactors, test-only
changes, and routine documentation cleanup unless project conventions require them.

Append one concise entry per logical change to the bottom of this section.
Use the fields below where useful; omit inapplicable fields:

### [vX.X.X] | YYYY-MM-DD
- **Type:** (feat | fix | change)
- **Description:** Resulting behavior and user impact.
- **Motivation:** Why was this change necessary?
- **Reference:** Optional issue, commit, or detailed project documentation.
- **Test Strategy:** How was this verified?

---

### [Unreleased] | 2026-11-20
- **Type:** feat
- **Description:** Add a callable, validated Fe2O3 run API and explicit versioned results with source/runtime provenance. Existing output is protected by default; legacy restarts/plots remain supported. Accepted short-case regression and maintained CPU/GPU/interface checks preserve existing physics and defaults. Saved MAT files now contain only `fe2o3_result` rather than the whole workspace.
- **Motivation:** Make runs reproducible and inspectable without changing numerical kernels or accepted trajectories.
- **Reference:** [Completion evidence and limitations](docs/cq01-05-verification.md).
- **Test Strategy:** Six migrated cases match accepted inputs and trajectories exactly; CPU/GPU and interface/restart/save checks pass, along with 23 expected corruption/path rejections.

---

### [Unreleased] | 2026-11-20
- **Type:** fix
- **Description:** Production saves now reject the acceptance manifest and its aliases even with explicit overwrite. Production and benchmark paths resolve relative to MATLAB's current directory after `cd`.
- **Motivation:** The manifest sat outside the protected baseline subtree, and Java's startup-directory resolution could misinterpret relative destinations.
- **Reference:** [CQ-04 follow-up verification](docs/cq04-followup.md).
- **Test Strategy:** Disposable manifest/path/link checks, tiny GPU interface/save/restart checks, and all 23 benchmark corruption/path checks pass; accepted artifacts and numerical kernels are preserved.

---

## Section 2: Error & Solution Log
For significant resolved errors, add a concise note here when useful. Skip
routine transient errors and avoid duplicating Section 1. Keep lengthy
investigations in issues or project documentation and link them here.

### [ERROR-ID] | Short Descriptive Title
- **Context:** Where/when did it happen?
- **Symptoms:** Brief description of the failure.
- **Root Cause:** Brief explanation of the cause.
- **Resolution:** Summary of the fix.
- **Prevention:** How to avoid recurrence, when useful.

---

## Section 3: Workflow for Agents
1. **Analyze** the task or error.
2. **Execute** code changes.
3. **Verify** with tests.
4. Update `CHANGELOG.md` for meaningful behavior changes or significant fixes,
   following the global changelog policy. Group related edits into one entry.
   Use Section 1 for changes and Section 2 for useful significant-error notes.
