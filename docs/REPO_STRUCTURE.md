# Repository Structure (Cleanup Baseline)

This document defines the current canonical layout to reduce path sprawl and make future refactors safer.

## Top-level conventions

- `app.R`: local development entrypoint.
- `inst/shiny/app.R`: package/runtime entrypoint used after install.
- `modules/`: actively maintained Shiny modules for repo-root app.
- `inst/shiny/modules/`: packaged mirror of modules used by installed app.
- `docs/`: design notes, bugfix history, and user/developer documentation.
- `scripts/`: operational scripts (install, debug, maintenance).
- `tests/`: automated tests and legacy tests.
- `archive/`: historical/legacy materials not in active runtime path.

## Active runtime paths

1. `app.R` + `modules/*`
2. `inst/shiny/app.R` + `inst/shiny/modules/*`

When changing one path, keep the other in sync unless deprecating intentionally.

## Script organization rules

- Keep executable utility scripts under `scripts/`.
- Keep one-off diagnostics in `scripts/debug_*` naming or `scripts/debug/` if added later.
- Avoid adding new ad-hoc scripts to repository root.

## Testing organization rules

- Keep automated Python checks in `tests/`.
- Keep R functional validations in `tests/` or clearly marked `tests/legacy/`.
- Prefer reproducible CLI-driven tests over manual-only notes.

## Cleanup plan (next phase)

- Consolidate root-level ad-hoc R check scripts into `scripts/` with stable names.
- Reduce duplicated documentation by linking from one index page.
- Add CI task to verify `modules/` and `inst/shiny/modules/` are synchronized.
