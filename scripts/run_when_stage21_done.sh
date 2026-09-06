#!/bin/bash
# scripts/run_when_stage21_done.sh
#
# Waits for the stage-21 K=4 fit to finish (posterior summary appears), then runs
# the K=4 residual & PCA analyses that append to the kcompare CSVs.
#
# Usage: bash scripts/run_when_stage21_done.sh
#        (safe to run in background: nohup bash scripts/run_when_stage21_done.sh > k4_analyses.log 2>&1 &)

set -e
cd "$(dirname "$0")/.."   # cd to project root

SUMMARY=outputs/tables/21_real_k4_t_posterior_summary.csv

echo "[watcher] waiting for $SUMMARY (poll every 60s)..."
until [[ -s "$SUMMARY" ]]; do sleep 60; done
echo "[watcher] $SUMMARY appeared. Running K=4 analyses..."

echo
echo "=== BFA_K=4 Rscript scripts/29b_kcompare_residuals.R ==="
BFA_K=4 Rscript scripts/29b_kcompare_residuals.R
echo
echo "=== BFA_K=4 Rscript scripts/31b_kcompare_pca.R ==="
BFA_K=4 Rscript scripts/31b_kcompare_pca.R

echo
echo "=== Cross-K summary tables after K=4 append ==="
echo "--- 29b_kcompare_summary.csv ---"
cat outputs/tables/29b_kcompare_summary.csv
echo
echo "--- 31b_kcompare_variance_shares.csv ---"
cat outputs/tables/31b_kcompare_variance_shares.csv
echo
echo "--- 31b_kcompare_reliable_counts.csv ---"
cat outputs/tables/31b_kcompare_reliable_counts.csv

echo
echo "[watcher] All K=4 analyses complete. Ready for PDF re-render."
