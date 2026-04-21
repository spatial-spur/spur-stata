---
name: Bug report
about: Report a reproducible problem with the spur Stata package
title: "[bug] "
labels: bug
---

## Summary

Describe the problem in 1-3 sentences.

## Minimal reproduction

Provide a copy-pasteable example if possible.

```stata

```

## Expected behavior

What did you expect to happen?

## Actual behavior

What happened instead? Include the full traceback, warning, or numerical mismatch if relevant.

```text

```

## Environment

- OS:
- Stata version:
- install method: `net install spur` / local ado path / manual copy into ado directory
- spur-stata version or commit:
- required dependency installed: `moremata` yes / no
- optional example dependencies installed, if relevant: `scpc` / `texdoc` / neither
- affected command or workflow: `spurtest` / `spurtransform` / `spurhalflife` / `example.do` / docs / installation

## Additional context

Anything else that might help reproduce or explain the issue.

- Does this reproduce on synthetic data, the Chetty example data, or only your own data?
- Does this affect `latlong` coordinates, Euclidean `s_*` coordinates, or both?
- If relevant, which mode or option is affected?
- Is the issue a crash, wrong result, parity mismatch, docs mismatch, install problem, or performance regression?
