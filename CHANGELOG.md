# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2025-12-25
### Added
- AST-only file-mode fallback: when dynamic import attempts fail (for local scripts with relative imports), LibInspector now performs a pure AST analysis and emits a Markdown/HTML file-level report.
- Multi-candidate (max-coverage) outputs: by default LibInspector attempts all import candidates (resolved package name, CLI entry module, constructed module name from file path, and raw file path) and generates a separate report per candidate.
- New API parameter: `run_all_candidates: bool` (default: `True`) to control multi-candidate behavior when calling `inspect_library()` programmatically.
- CLI flag: `--no-multi-candidate` / `--multi-candidate` to control multi-candidate runs from the command line.
- Report header enhancements: each generated report now includes discovery context (`discovered_via`), tried-candidates table (with errors), used-candidate hint, package version (best-effort), run-mode note, and a short `sys.path` snapshot to aid debugging.

### Changed
- CLI entry-point candidate ordering adjusted: CLI entry modules are appended behind package/module candidates so single-run mode prefers package/module imports.
- CLI behavior: the CLI currently triggers multi-candidate runs by default (per-candidate Markdown/HTML files are created). If you prefer a single consolidated run programmatically, call `inspect_library(..., run_all_candidates=False)`.
- Output filenames: when multiple candidates are produced, output files are named by appending `_candidate_<n>_<safe-name>` to your requested output path (e.g. `report_candidate_1_TensorVis.cli.md`).

### Fixed
- Improved import candidate resolution and added safer sys.path injection heuristics.



## [0.1.0] - 2025-11
### Added
- Initial release supporting dynamic inspection of installed packages and reverse lookup of CLI commands (resolve CLI -> package -> entry module), producing Markdown + HTML reports with dependency graphs and function flow visualizations.

---
Notes:
- For full details and examples see the README and tests/ outputs in `/test/`.
