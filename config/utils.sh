#!/bin/bash
# ====================================================
# 📦 Shared Utility Functions for CrispAstro-Seq Pipeline
# ====================================================
  
# 🕒 Initialize general timestamps (used for logs and output naming)
initialize_timestamps() {
  TIMESTAMP=$(TZ=UTC date +%y%m%d_%H%M)
  START_TIME_KSA=$(TZ=Asia/Riyadh date +'%Y-%m-%d %H:%M AST (%:z)')
  START_TIME_UTC=$(TZ=UTC date +'%Y-%m-%d %H:%M UTC (%:z)')
}

# 🔧 Check tool presence silently and set TOOL_VERSION
# Usage: check_tool_version TOOL_NAME
check_tool_version() {
  local tool_name="$1"

  if ! command -v "$tool_name" &>/dev/null; then
    echo "❌ ERROR: Required tool '$tool_name' is not installed or not in PATH."
    exit 1
  fi

  TOOL_VERSION=$("$tool_name" --version 2>&1 | head -n 1)
}

# ⏱ Capture runtime and end time based on START_TIMESTAMP (set in main script)
# 🧾 Final runtime report (end time + duration) — safe one-call version
report_runtime() {
  if [[ -z "${START_TIMESTAMP:-}" ]]; then
    echo "❌ ERROR: START_TIMESTAMP is not set. Set it before calling report_runtime()."
    exit 1
  fi

  END_TIMESTAMP=$(date +%s)
  RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
  RUNTIME_FMT=$(date -ud "@$RUNTIME" +'%H hrs %M min %S sec')
  END_TIME_KSA=$(TZ=Asia/Riyadh date +'%Y-%m-%d %H:%M:%S AST (%:z)')
  END_TIME_UTC=$(TZ=UTC date +'%Y-%m-%d %H:%M:%S UTC (%:z)')

  echo ""
  echo "✅ Finished at       : $END_TIME_KSA | UTC: $END_TIME_UTC"
  echo "⏱️ Total Runtime     : $RUNTIME_FMT"
}

