#!/bin/bash
# Periodically sync Modal Volume results to local machine.
# Usage: ./scripts/sync_modal_results.sh [interval_seconds]
#   Default interval: 900 (15 minutes)

VOLUME="rfab-cancer-drivers-results"
LOCAL_DIR="./results/cancer_drivers_modal_sync"
INTERVAL="${1:-900}"

mkdir -p "$LOCAL_DIR"

echo "Syncing Modal Volume '$VOLUME' → $LOCAL_DIR every ${INTERVAL}s"
echo "Press Ctrl+C to stop"
echo ""

while true; do
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[$ts] Syncing..."

    # Sync batch state
    modal volume get "$VOLUME" .batch_state.json "$LOCAL_DIR/.batch_state.json" --force 2>/dev/null
    if [ -f "$LOCAL_DIR/.batch_state.json" ]; then
        echo "  batch_state.json synced"
        python3 -c "
import json, sys
state = json.load(open('$LOCAL_DIR/.batch_state.json'))
campaigns = state.get('campaigns', {})
for name in sorted(campaigns):
    info = campaigns[name]
    status = info.get('status', '?')
    stage = info.get('current_stage', '')
    elapsed = info.get('elapsed_seconds', 0)
    print(f'  {name:<30} {status:<12} {stage:<25} {elapsed/3600:.1f}h')
print(f'  Total: {len(campaigns)} campaigns tracked')
" 2>/dev/null
    fi

    # Sync logs for all campaigns
    for campaign_dir in $(modal volume ls "$VOLUME" cancer_drivers/ 2>/dev/null | grep "cancer_drivers/"); do
        campaign=$(basename "$campaign_dir")
        mkdir -p "$LOCAL_DIR/$campaign/logs"
        modal volume get "$VOLUME" "cancer_drivers/$campaign/logs/" "$LOCAL_DIR/$campaign/logs/" --force 2>/dev/null
    done
    echo "  Logs synced"

    # Sync batch summary if exists
    modal volume get "$VOLUME" batch_summary.json "$LOCAL_DIR/batch_summary.json" --force 2>/dev/null

    echo "[$ts] Sync complete. Next sync in ${INTERVAL}s"
    echo ""
    sleep "$INTERVAL"
done
