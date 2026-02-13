#!/bin/bash
# Quick progress check for all 10 cancer driver campaigns on Modal.
# Usage: ./scripts/check_progress.sh

VOLUME="rfab-cancer-drivers-results"
CAMPAIGNS="b7h3_vhh cd47_vhh ceacam5_vhh egfr_cetuximab_vhh egfrviii_vhh epha2_vhh gpc2_vhh her2_domIV_vhh msln_cterm_vhh msln_nterm_vhh"
TARGET=500
TMPDIR="${TMPDIR:-/tmp}/rfab_progress_$$"
mkdir -p "$TMPDIR"

echo ""
echo "As of: $(date '+%Y-%m-%d %H:%M:%S %Z') ($(date -u '+%H:%M UTC'))"
echo ""
printf "%-26s %8s %8s %8s\n" "Campaign" "Designs" "% Done" "Stage"
echo "------------------------------------------------------"

total_designs=0
for campaign in $CAMPAIGNS; do
    logfile="$TMPDIR/${campaign}.log"
    modal volume get "$VOLUME" "cancer_drivers/${campaign}_v1/logs/${campaign}.log" "$logfile" --force >/dev/null 2>&1

    if [ -f "$logfile" ] && [ -s "$logfile" ]; then
        designs=$(grep -c "Finished design" "$logfile" 2>/dev/null || echo 0)
        pct=$((designs * 100 / TARGET))

        if grep -q "Stage 3" "$logfile" 2>/dev/null; then
            stage="RF2"
        elif grep -q "Stage 2" "$logfile" 2>/dev/null; then
            stage="ProteinMPNN"
        elif grep -q "Stage 1" "$logfile" 2>/dev/null; then
            stage="RFdiffusion"
        else
            stage="starting"
        fi

        # Check if campaign completed
        if grep -q "Checkpoint saved: stage3" "$logfile" 2>/dev/null; then
            stage="DONE"
        fi

        total_designs=$((total_designs + designs))
        printf "%-26s %8d %7d%% %8s\n" "$campaign" "$designs" "$pct" "$stage"
    else
        printf "%-26s %8s %8s %8s\n" "$campaign" "-" "-" "no log"
    fi
done

echo "------------------------------------------------------"
printf "%-26s %8d\n" "Total designs" "$total_designs"
echo ""

# Runner stats
python3 -c "
import modal
f = modal.Function.from_name('rfantibody-harness-cancer-targets', 'run_campaign')
s = f.get_current_stats()
print(f'GPU runners: {s.num_total_runners} active, backlog={s.backlog}')
" 2>/dev/null

rm -rf "$TMPDIR"
