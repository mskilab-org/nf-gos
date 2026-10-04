# JaBbA / allelic-balance container optimization handoff

Date: 2026-10-04. Repository on the originating machine: `/gpfs/home/hadik01/git/nf-gos`.

## Goal and current status

Move the runtime formulation patches into maintained gGnome/JaBbA package source, build a reproducible replacement Docker image, verify the actual solver path, then remove the obsolete runtime patches in one coordinated pipeline cutover. This document is the handoff; no Docker image was built or published during this change.

The pipeline currently works without rebuilding its images:

| Process | Declared image | Active implementation |
| --- | --- | --- |
| `JABBA` | `mskilab/jabba:0.0.9` | Staged `bin/jabba_optimized.R` patches the `balance` import in JaBbA's imports environment, then sources the packaged `extdata/jba` launcher. `jbaLP` resolves that patched import. |
| `NON_INTEGER_BALANCE` | `mskilab/jabba:0.0.8` | New staged `bin/non_integer_balance_optimized.R` sources staged `bin/non_integer_balance.R` in a private environment whose `balance` function uses a patched copy of gGnome's implementation. |
| `LP_PHASED_BALANCE` | `mskilab/jabba:0.0.8` | `bin/lp_phased_balance.R` defines its own large `balance()` implementation. Updating the installed gGnome package alone does not replace it. |

Neither runtime wrapper changes the installed gGnome namespace or persists patches across R processes. Merely switching the non-integer image to `0.0.9` does not inherit the JaBbA wrapper's patch.

## Changes made in this handoff's originating session

- Added `bin/non_integer_balance_optimized.R`, following the existing JaBbA wrapper's guarded deparse/substitution approach.
- Added `--threads` (standalone default 1) and `--mipemphasis` (default 0) to the algorithm's CLI. These options are consumed by the optimized wrapper; invoking the base script directly does not apply the optimization.
- Updated `modules/local/allelic_cn/main.nf` to stage both R files, invoke the wrapper, pass `${task.cpus}`, and pass emphasis. Script content now participates in Nextflow task caching rather than being found only through `NEXTFLOW_PROJECT_DIR`.
- Updated `subworkflows/local/allelic_cn/main.nf` with the two script channels and emphasis argument. Its public workflow inputs are unchanged.
- Added `params.mipemphasis_non_integer_balance = 0` in `nextflow.config`. Use `--mipemphasis_non_integer_balance 1` for a feasibility-focused run. Default gap/time limit remain 0.001/3600 seconds; thread count follows the effective process allocation, not a hardcoded 16.
- Corrected the affected process's stub filename to `non_integer.balanced.gg.rds`, matching its existing output declaration.
- Images, JaBbA wrapper behavior, phased formulation, biological filtering, objective weights, and existing CN >999 masking were not changed.

## Exact non-integer runtime behavior to migrate

### Formulation

The patched copy of `gGnome::balance()`:

1. After `gg = gg$copy`, computes `observed.cn = max(abs(gg$nodes$dt$cn), na.rm=TRUE)`, using 0 if the result is nonfinite; sets `M = max(50, ceiling(2 * observed.cn))` and logs it.
2. Immediately after the existing loose-end bound initialization, sets loose-in/out upper bounds to M and all `vtype == "B"` bounds to `[0,1]`.
3. Leaves the subsequent `nonintegral` block untouched. It restores terminal-slack lower bounds to `-0.4999`, gives REF edges that same lower bound, and retains chromosome slush bounds `[-0.49999, 0.49999]`.
4. Adds an R-level `mipemphasis` argument and sets `control$mipemphasis` immediately before the `Rcplex2` call.

The script is still a MILP: `lp=TRUE` is the linear objective mode, not continuous relaxation. Node and ALT-edge CN remain integer; REF-edge/terminal CN and chromosome offsets can be continuous. Do not replace all variables with continuous variables or impose blanket nonnegative bounds after the nonintegral block.

The existing JaBbA wrapper uses the same observed-CN M heuristic and binary/loose-end bounds, but patches JaBbA's import instead of the non-integer script's private environment. It does not add emphasis or thread controls.

**M is a heuristic cap, not proven formulation equivalence.** It also bounds delta variables, not only big-M indicator coefficients. Assess high-CN samples, fixed nodes/edges, unknown CN, and forced junctions before generalizing this rule into a package-wide default. Respect explicit M supplied by other callers in the maintained API; the current wrappers unconditionally replace M. Do not silently erase observations or add narrow per-node MLE windows to force feasibility.

Explicit residual bounds were not added by the new wrapper. In this LP path the delta constraints already bound constrained residuals. Do not blindly copy the phased script's residual/slack rules into the non-integer path.

### Solver settings and precedence

- Verified solver runtime: CPLEX **22.1.2.0** in the local `0.0.8` environment with host CPLEX mounted at `/opt/cplex`. gGnome reports version `0.1`, which is not enough to identify source provenance.
- The installed `Rcplex2` does not support thread count in its control list. The wrapper writes `non_integer.cplex.prm`, sets `ILOG_CPLEX_PARAMETER_FILE`, and lets CPLEX read it at `CPXopenCPLEX`.
- Existing parameter-file settings are retained except old/modern thread and emphasis entries, which are removed. The requested thread value is appended. Without an inherited file the header is `CPLEX Parameter File Version 22.1.2.0`, matching the verified deployment. This assumption must be revisited for another runtime.
- **Emphasis must go through `control$mipemphasis`.** An actual solve showed that specifying emphasis only in the parameter file was reset to balanced mode by Rcplex2's subsequent control processing. Threads remained honored. This was corrected before delivery.
- The wrapper validates positive integral thread counts and emphasis in 0..4 for CPLEX. Balanced emphasis 0 remains the default. Warm starts and external restarts were not introduced.
- The optional Gurobi branch retains the formulation patch but does not use these CPLEX settings. Its installed `run_gurobi()` accepts a `threads` formal but does not forward it to Gurobi parameters. Gurobi was not exercised; do not claim its resource control is fixed by this change.
- Do not infer the linked CPLEX version from a `cplex` executable on PATH. Earlier observations had CLI 12.8 versus linked runtime 22.1.2.

## Work for the container-build agent

1. **Locate and record actual build provenance.** This nf-gos checkout has no Dockerfile/Containerfile/Singularity recipe. Locate the maintained image-build repository on your machine or inspect image labels. Record image digest, R/Bioconductor versions, package source commits, and linked solver library/version. Do not invent a Dockerfile path or rebuild from floating package branches.
2. **Patch maintained gGnome source.** Find the implementation of `balance`, update bound construction at the source, preserve nonintegral semantics, and expose the necessary solver controls. Choose an explicit adaptive-M policy/API that does not override unrelated callers' deliberate bounds. Patch compiled Rcplex parameter handling only if moving thread control out of the parameter-file route; that route already works without a native change.
3. **Audit JaBbA integration.** Confirm `ramip_stub -> jbaLP -> balance -> Rcplex2` in the selected package version. JaBbA's imported function must resolve the updated implementation. Audit its own Rcplex2 copy separately if you change native parameter handling; JaBbA and gGnome expose different bindings.
4. **Treat LP phased balance separately.** Compare its script-local implementation with the packaged version before any replacement. It contains extra phasing/CNLOH logic; deleting the local function merely because gGnome was upgraded is not a valid migration.
5. **Build an immutable new linux/amd64 image.** Pin source commits/base digest. Follow the lab's CPLEX licensing/distribution policy. The repository documents CPLEX as proprietary and separately installed; this deployment mounts host CPLEX. Do not publish solver binaries, license files, or credentials without authorization. Do not overwrite the existing 0.0.8/0.0.9 tags.
6. **Verify the new image before switching the pipeline.** Use the acceptance checklist below and archive solver logs, package provenance, objective/bound/elapsed time, and output validation.
7. **Perform a clean cutover.** Replace wrapper-supplied formulation behavior with calls to the maintained package API, remove obsolete patch functions/import mutation, migrate every affected caller and staged input/channel, and retain script staging for the algorithm files still used. Preserve effective thread/emphasis behavior. Update both Docker/Singularity image branches in affected modules. Do not indiscriminately upgrade unrelated processes in the same files. Update this document/changelog with the final image tag/digest and source commits.

## Acceptance checklist

- Run actual JaBbA LP/MILP and non-integer balance paths, not just package imports or `--help`.
- Verify generated model domains: binary indicators `[0,1]`; node and ALT CN integral; REF/terminal fractional lower bounds preserved; chromosome slush `[-0.49999,0.49999]`; loose-end upper bounds consistent with the selected M policy.
- Confirm CPLEX logs match allocated threads at two different allocations and that emphasis 0 and 1 actually select balanced/feasibility modes. Include an inherited parameter file with conflicting settings to check precedence.
- Confirm feasible solutions satisfy the actual sparse constraint matrix and variable bounds numerically, and validate output graph balance/integrality as appropriate. Account for the existing rounded graph outputs and separately stored slush metadata.
- Cover default and high-CN graphs, forced ALT edges, fixed bounds, missing CN, chromosome tips, and nonstandard contigs. The earlier raw-graph smoke with nonstandard contigs hit `NA's in (i,j) are not allowed` before the solve; do not claim this was repaired. Investigate it independently rather than silently dropping contigs in production.
- With equal time/thread/seed settings and production coverage-derived weights, compare baseline, bounds-only, smaller-M, and emphasis variants separately. Report objective and best bound, not relative gap alone. A better incumbent can coexist with a wider gap.
- Confirm full Nextflow process outputs and resume/cache behavior after changing staged scripts. A stub run is only a channel/filename check, not solver proof.
- Fail clearly on an unsupported package implementation while runtime patches remain. After the package cutover, remove obsolete runtime patching rather than applying it twice or retaining a silent fallback.

## Verification already performed here

All smoke inputs and outputs were isolated from production; temporary scripts were removed after verification. No production graph was overwritten.

- Current `jabba_optimized.R` patch setup in image 0.0.9: `jbaLP` resolved the patched import; installed gGnome namespace remained unchanged.
- New wrapper `--help` executed in image 0.0.8.
- Instrumented **real** CPLEX solve using the available WG-26-125 standard-chromosome subgraph, unit node weights, unknown edge CN, ALT lb=1: 22,438 variables / 21,731 rows. Checked all relevant variable domains, M=50, emphasis 1, two threads, and retained inherited random seed 17. Feasible incumbent objective 515.0001 at five seconds; maximum constraint violation approximately `1.35e-8`. Checked unsupported-formulation and invalid-parameter failures. This is a behavior smoke, not a production performance claim.
- Full optimized CLI ran through the original algorithm using a prepared binstats cache, tiny indexed reference, `allin=FALSE`, no hets, emphasis 0, two threads. It completed balance and downstream conversion and wrote both graph RDS outputs.
- Nextflow 25.10.4: stub and **actual containerized `NON_INTEGER_BALANCE` process** both completed and emitted `non_integer.balanced.gg.rds` and `hets.gg.rds`. Actual task logged M=50, two CPLEX threads, balanced emphasis, and `CPXMIP_TIME_LIM_FEAS`. Output objective was 463.0001, relative gap 0.393497 at five seconds. The real module smoke used the same prepared binstats cache and tiny reference; coverage fitting, junction mappability filtering, and het inference were not exercised.
- Existing package startup warnings about Bioconductor imports and expired Ensembl SSL certificates were observed; they did not prevent these runs. They were not suppressed as a fix.

## Prior findings to retain

The LP phased-balance controlled comparison on WG-26-107 improved the incumbent from 91,404 to 67,928 in the same 3,600 seconds, driven by tighter M plus feasibility emphasis. Warm-start comparison did not improve the incumbent; default remains off. Narrow per-node CN windows were rejected after infeasibility. Do not change the known CNLOH reward-coefficient precedence or other dead flags as part of this container migration without a separate decision.

Correction to an earlier explanation: malformed binary bounds do not establish that CPLEX treated every B variable as a general integer. The baseline non-integer runtime still reported thousands of binaries, and its first presolved binary/general counts matched the explicitly bounded model. Explicit bounds removed the rounding warning; the warning alone is not a performance mechanism or benchmark.
