# Environments

Two conda environments, both built against R 4.5.

| File | Environment | Use it for |
|---|---|---|
| `matrix2-core.yml` | `matrix2-core` | The method itself: `shared_edges.r`, `simple_code.r`, `test_edge_conservation.r`. Small, solves quickly. |
| `matrix2-expression.yml` | `matrix2-expression` | Real data: transcript import, differential expression, GO enrichment. A superset of core. |

```bash
mamba env create -f envs/matrix2-core.yml
conda activate matrix2-core
Rscript test_edge_conservation.r
```

`test_edge_conservation.r` needs only the `Matrix` package, so it also runs under a bare R installation without either environment.

## Checking an environment works

```bash
conda activate matrix2-expression
Rscript envs/smoke_test.R
```

Loading a package is not the same as it working — a bioconductor package built against the wrong R version imports fine and then fails inside compiled code. `smoke_test.R` runs DESeq2, apeglm shrinkage, edgeR, topGO and tximport on tiny synthetic data with a planted effect, so a broken build surfaces immediately rather than three hours into an analysis.

## Why R is pinned but the packages are not

Every `bioconductor-*` package is built against one specific R minor version — the `r45` in a build string like `r45hdfd78af_0`. Mixing packages built for different R versions does not resolve, so `r-base=4.5` is pinned in both files. The individual package versions are left open so the environments stay installable as bioconda rebuilds them.

That is the right trade-off for a file you install from, but it is **not** reproducibility. For that, export a lock once the analysis is settled:

```bash
conda activate matrix2-expression
conda env export --no-builds > envs/matrix2-expression.lock.yml
```

Commit the lock alongside the results it produced, and record `sessionInfo()` with each set of figures. Reviewers ask which DESeq2 version produced a result, and "the one bioconda served that day" is not an answer.

## Adding a package

Add it to the `.yml` rather than installing into a live environment, then recreate:

```bash
mamba env remove -n matrix2-expression
mamba env create -f envs/matrix2-expression.yml
```

An environment that has been hand-patched with `install.packages()` cannot be rebuilt by anyone else, including you in six months.
