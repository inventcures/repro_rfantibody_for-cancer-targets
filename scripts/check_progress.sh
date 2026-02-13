#!/bin/bash
# Quick progress check for all 10 cancer driver campaigns on Modal.
# Tracks all 3 pipeline stages: RFdiffusion → ProteinMPNN → RF2
# Shows RF2 (Stage 3) prediction progress as scored/total.
# Usage: ./scripts/check_progress.sh

VOLUME="rfab-cancer-drivers-results"
CAMPAIGNS="b7h3_vhh cd47_vhh ceacam5_vhh egfr_cetuximab_vhh egfrviii_vhh epha2_vhh gpc2_vhh her2_domIV_vhh msln_cterm_vhh msln_nterm_vhh"
TMPDIR="${TMPDIR:-/tmp}/rfab_progress_$$"
mkdir -p "$TMPDIR"

echo ""
echo "As of: $(date '+%Y-%m-%d %H:%M:%S %Z') ($(date -u '+%H:%M UTC'))"
echo ""
printf "%-24s %6s %6s %12s %14s\n" "Campaign" "S1" "S2" "S3 (RF2)" "Stage"
echo "--------------------------------------------------------------------"

done_count=0
total_scored=0
total_expected=0

for campaign in $CAMPAIGNS; do
    dir="cancer_drivers/${campaign}_v1/pipeline"
    logfile="$TMPDIR/${campaign}.log"

    s1="-"; s2="-"; s3_display="-"; stage="waiting"
    scored=0; expected=0

    # Stage detection via volume ls (fast, no downloads)
    files=$(modal volume ls "$VOLUME" "$dir/" 2>/dev/null)

    if echo "$files" | grep -q "01_backbones.qv" 2>/dev/null; then
        s1="yes"
        stage="RFdiffusion"
    fi
    if echo "$files" | grep -q ".checkpoint_stage1" 2>/dev/null; then
        s1="done"
    fi
    if echo "$files" | grep -q "02_sequences.qv" 2>/dev/null; then
        s2="yes"
        stage="ProteinMPNN"
    fi
    if echo "$files" | grep -q ".checkpoint_stage2" 2>/dev/null; then
        s2="done"
    fi

    has_predictions=false
    if echo "$files" | grep -q "03_predictions.qv" 2>/dev/null; then
        has_predictions=true
        stage="RF2"
    fi
    if echo "$files" | grep -q ".checkpoint_stage3" 2>/dev/null; then
        stage="DONE"
        done_count=$((done_count + 1))
    fi

    # RF2 progress: download log and count completed predictions
    if [ "$has_predictions" = true ] || [ "$stage" = "DONE" ]; then
        modal volume get "$VOLUME" "cancer_drivers/${campaign}_v1/logs/${campaign}.log" "$logfile" --force >/dev/null 2>&1
        if [ -f "$logfile" ] && [ -s "$logfile" ]; then
            scored=$(grep -c "\[RF2\] Completed:" "$logfile" 2>/dev/null || echo 0)
            # Count ProteinMPNN successes to get total sequences (= RF2 expected)
            mpnn_count=$(grep -c "reported success" "$logfile" 2>/dev/null || echo 0)
            expected=$((mpnn_count * 5))
            if [ "$expected" -eq 0 ]; then
                expected="?"
                s3_display="${scored}/?"
            else
                pct=$((scored * 100 / expected))
                s3_display="${scored}/${expected}"
                if [ "$stage" = "DONE" ]; then
                    s3_display="done"
                fi
            fi
        fi
    fi

    total_scored=$((total_scored + scored))
    if [ "$expected" != "?" ] 2>/dev/null; then
        total_expected=$((total_expected + expected))
    fi

    printf "%-24s %6s %6s %12s %14s\n" "$campaign" "$s1" "$s2" "$s3_display" "$stage"
done

echo "--------------------------------------------------------------------"
if [ "$total_expected" -gt 0 ]; then
    total_pct=$((total_scored * 100 / total_expected))
    printf "%-24s %6s %6s %12s %14s\n" "TOTAL" "" "" "${total_scored}/${total_expected}" "${done_count}/10 done"
else
    printf "%-24s %6s %6s %12s %14s\n" "TOTAL" "" "" "${total_scored}" "${done_count}/10 done"
fi
echo ""

# Runner stats
python3 -c "
import modal
f = modal.Function.from_name('rfantibody-harness-cancer-targets', 'run_campaign')
s = f.get_current_stats()
print(f'GPU runners: {s.num_total_runners} active, backlog={s.backlog}')
" 2>/dev/null

rm -rf "$TMPDIR"
