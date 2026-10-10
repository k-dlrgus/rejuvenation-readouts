# CLAUDE.md — standing rules for this repository

This project's method is pre-registration. Its credibility rests on the record
being honest about what was decided when, so these rules are about the record
as much as about the code.

## Scope of a session

- **One task per session.** A task is one pre-registered analysis. Do not start
  a second analysis because the first finished early.
- **Two hand-offs only.** The pre-registration goes to chat before anything
  runs; results go to chat after. Nothing in between.
- **Hand over files, never summaries.** A summary in the chat is not a result.
  Write the findings file and the tables, then point at them.
- **Never resume a stopped run.** If a STOP fires, write the findings with the
  STOP key and exit. Do not compute after a STOP, and do not restart it with an
  adjustment.

## Pre-registration

- The flag file is written **before any number exists**, with a timestamp, and
  reproduced **verbatim** at the top of the findings file.
- It fixes, in advance: the statistic, the nulls, the positive controls, the
  decision thresholds, and **the exact wording of each possible conclusion**.
- Only the pre-registered reading whose conditions fired is reported as the
  conclusion. Diagnostics and descriptive summaries are not gated results and
  are labelled as such.
- **If the flag is ambiguous or self-contradictory, stop and ask.** Do not
  resolve it by picking the reading that lets the run proceed. If a definition
  is genuinely missing, fix it in writing before computing and say that you did.
- **Disclose every deviation**, including ones that look like housekeeping:
  where you wrote, what you imported, what you added beyond the flag, what you
  decided that the flag did not specify.
- A reading that is superseded is **marked superseded, never edited**. The old
  text stays.

## Before computing anything new

- **Reproduce the stored numbers first**, with this task's code, and report the
  comparison as a table. Counts exact; numerics within the tolerance the flag
  names. If reproduction fails, STOP and name every mismatch.
- Frozen instruments — the ruler, the anchors, any fitted axis — are **read, not
  refit**, unless the flag says otherwise. If something must be refit, say so
  prominently and state what entered the fit.
- Never write into a previous task's results directory. New task, new directory.

## Honesty obligations specific to this project

- **Report the chance baseline beside any count.** A fifth-percentile floor
  passes 5% of true nulls by construction; a bare count without that baseline is
  misleading.
- **Check any shortlist against cell count** before reporting it. Group size in
  this screen is an accident of the experiment, and a previous shortlist turned
  out to be the factors that happened to capture the most cells.
- **A null is only a null if the statistic could have seen the effect.** Report
  the recovery power of the statistic at the depth in question alongside any
  zero.
- **Circularity is disclosed, not refitted away.** If an axis, basis or
  reference was estimated from the same cells it is later applied to, say so
  where the result appears.
- State what a null does **not** cover. Results on neonatal or passaged cells
  with no aged donor do not transfer to aged cells, and the findings should say
  so.

## Repository hygiene

- Layout per task: `src/<task>_run.py`, `results/<task>/`, `FINDINGS_<TASK>.md`.
- **New commits, never amend.** One commit per task, in dependency order, each
  carrying its findings file, its `src/` entry and its results directory. The
  committed code must run and its gates must be re-verifiable from tracked files
  alone.
- Git-ignore regenerable binaries of 1 MB or more; they are deterministic from
  seeds. Keep result CSVs tracked. Note in a comment where an ignored raw input
  came from.
- **Do not push.** The repository is public and publication timing is a
  strategy decision made outside the session. Commit locally and stop.
- Do not commit this file's sibling raw data under `data/raw/`.

## Tone

Report what happened, including the parts that undercut the result. The session
that did the work is the worst placed to judge it, so make the findings file
good enough that someone else can.
