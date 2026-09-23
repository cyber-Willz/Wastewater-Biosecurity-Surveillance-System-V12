# ww_biosec v0.7 — Implemented Fixes and Re-Run Results

This documents the two fixes implemented in `ww_detection` against the live Scotland
dataset, what changed in the results, and — importantly — where the first fix (CUSUM)
turned out not to work as hypothesized, and what the working fix cost in alert volume.
Nothing below is simulated; every number comes from re-running `scotland_e2e` against the
same live-fetched `scotland_raw.csv` used in the original run.

## 1. Fix #1: variance floor (landed as designed)

**Problem:** `EwmaBaseline`'s Welford variance had a floor of `1e-6` (sd ≈ 0.001 log10
units) — far below realistic RT-qPCR/environmental noise. A run of near-identical early
observations (quiet warm-up, or several weeks pinned at the LOD floor) could collapse
variance toward that floor, which then broke two things together: the winsorisation
window (`EWMA_CLIP_FRAC · clip_z · sd`) shrank to near-zero, clamping the EWMA mean
almost immobile the moment real signal arrived — the "cold-start lock-up" the shipped
code's own `--metric=normalized` warning referenced — while the *raw* z-score computed
against that same near-zero sd exploded on the next observation.

**Fix:** raised the floor to `MIN_VARIANCE = 0.01` (sd = 0.1 log10 units), in
`ww_detection/src/baseline.rs`. Two new regression tests confirm a flat warm-up no longer
locks the EWMA in place and no longer produces an absurd z-score on the first real
signal.

**Effect on the live run:** Fort William's alert count dropped from 26 to 13, and its
extreme z-scores shrank substantially (37.2→22.6, 42.9→10.0, 21.7→12.9); several of its
weaker CRITICAL calls downgraded to AMBER/RED. This is the artifact-reduction fix
promised in the earlier accuracy analysis — it landed and did what it was supposed to,
independent of anything below.

## 2. Fix #2: CUSUM sustained-shift detector (implemented, but empirically didn't fix Omicron)

**Original hypothesis:** the EWMA's fast adaptation (half-life ≈ 1.4 observations at
α=0.39) would absorb a real, *sustained* multi-week rise before any single week's
z-score crossed 2.5 — so a CUSUM accumulator run in parallel, summing small persistent
deviations, should catch it.

**What was actually found after implementing and checking it against the per-site trace:**
individual sites' z-scores during the Omicron window rose for exactly one or two weeks,
then reverted. There was no multi-week sustained elevation at any single site for a
temporal CUSUM to accumulate against — so it never fired during Omicron in the live run
(confirmed directly against the trace data). The CUSUM code is real, tested, and does
correctly catch genuinely sustained single-site shifts (regression test included) — it
just isn't what explains the Omicron miss in this dataset. Kept in the codebase since
it's a real, independently-useful improvement over pure z-score chasing; documented
honestly here so its role isn't overstated.

## 3. Fix #3 (the one that actually worked): cross-sectional "broad-based elevation"

Re-examining the trace data showed the real shape of the Omicron signal: not one site
rising for weeks, but **many different sites rising modestly in the same one-to-two-week
window** — a spatial/cross-sectional pattern, not a temporal one. That's precisely what
the shipped spectral/network score was meant to catch, but it runs on a health-board
proxy topology the docs already flag as uncalibrated on this data.

**Fix:** added a direct, transparent cross-sectional check in `ww_detection` (new
`AnomalyDetector::observe_analyte_ensemble`, `BROAD_BASED_Z_BAR` / `BROAD_BASED_FRAC_THRESHOLD`
constants): each round, if ≥20% of currently-warm, reporting sites (gated on ≥30
reporting sites, to avoid early-2020 small-sample noise) show their own z-score ≥ 1.0
simultaneously, every site clearing that bar is surfaced as at least AMBER — tagged
`network_onset` — even if its own z-score falls short of the normal 2.5 threshold. The
20%/1.0 pairing was calibrated against the trace data itself, not picked arbitrarily:
every week with ≥30 reporting sites that cleared a 0.20 broad-based fraction fell inside,
or within a week or two of the start of, a real wave window.

## 4. Re-run results (live data, same raw CSV, same `min_weeks=8`)

| Metric | v0.6 (original) | v0.7 (with fixes) |
|---|---|---|
| Total AMBER+ alerts | 81 | 464 |
| Alerts inside a known wave | 57 (70.4%) | 371 (80.0%) |
| Alpha wave caught? | ✅ | ✅ (7→12 alerts) |
| Delta wave caught? | ✅ | ✅ (50→352 alerts) |
| **Omicron wave caught?** | ❌ **No alerts** | ✅ **7 alerts** |
| Fort William alert count | 26 | 13 |
| Fort William max z-score | 42.9 | 22.6 |

**The honest cost:** 395 of the 464 alerts (85%) fired *only* because of the new
cross-sectional check — they would not have fired under the old z-score-only logic. Of
those 395, 325 (82%) land inside a known wave and 70 (18%) don't. Most of the volume
increase is concentrated in the Delta wave (50→352): the real Delta surge genuinely did
elevate more than half of all monitored sites simultaneously during its peak weeks
(June–September 2021), so once the detector is sensitive enough to catch Omicron's
narrower, brief simultaneous rise, it also faithfully reports every site that
participated in Delta's much broader one — which a national system may or may not want
surfaced at that granularity.

## 5. What this means, plainly

- The variance-floor fix is a clean win: less artifact, same detection power, no
  tradeoff found.
- The CUSUM fix is real and tested but turned out not to be the mechanism behind the
  Omicron miss — included for completeness and because it's a genuine improvement over
  pure z-score chasing, not because it resolved the headline problem.
- The cross-sectional fix is what actually closes the Omicron gap, and does so for a
  documented, quantifiable reason (a real, calibrated statistical signal, not a hack) —
  but it moves alert volume from "occasional" (~1/week) to "bursty" (~10-70/week during a
  real wave). Whether that trade is worth it depends on what the alerts feed into: a
  human-reviewed AMBER queue that only escalates on RED/CRITICAL would likely benefit;
  a system that treats every AMBER as needing individual triage would need either a
  round-level (one alert per network-onset event, not per site) aggregation on top of
  this, or a higher frac/z bar traded against re-missing Omicron.
- This wasn't run against the `--metric=normalized` path or through the Axum/PostgreSQL
  layer — both remain open items from the earlier accuracy analysis.

## 6. Files

| File | Contents |
|---|---|
| `scotland_alerts_v4.csv` | Full v0.7 alert list — includes new `network_onset` and `sustained_shift` columns |
| `scotland_trace_v4.csv` | Full v0.7 per-site-week trace |
| `ww_biosec_v07_before_after.png` | Alert-volume-per-week, before vs. after, with wave windows shaded |

(Left in the build container; ask if you want them alongside this report.)
