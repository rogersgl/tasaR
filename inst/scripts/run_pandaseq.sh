#!/usr/bin/env bash
# Bash version for Linux
# Place in same folder as Input.csv and unpaired reads
# DO NOT CHANGE COLUMN ORDER in Input.csv (is crude and checks positions, not column names)

set -euo pipefail

# Get script directory (portable)
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_FILE="${DIR}/Input.csv"
LOGFILE="${DIR}/run_pandaseq.log"

# Start fresh log
: > "$LOGFILE"
echo "Start: $(date)" >> "$LOGFILE"

mkdir -p "$DIR/logs" "$DIR/merged"

# Binary to call (use full path if needed)
PANDASEQ_BIN="pandaseq"
if ! command -v "$PANDASEQ_BIN" >/dev/null 2>&1; then
  # Try a common location; comment/uncomment as needed:
  if [ -x "/usr/local/bin/pandaseq" ]; then
    PANDASEQ_BIN="/usr/local/bin/pandaseq"
  fi
fi
echo "Using pandaseq: $(command -v "$PANDASEQ_BIN" 2>/dev/null || echo 'NOT FOUND')" >> "$LOGFILE"

# Ensure input exists
if [ ! -f "$INPUT_FILE" ]; then
  echo "ERROR: Input file not found: $INPUT_FILE" | tee -a "$LOGFILE"
  read -n1 -r -p $'\nPress any key to close...\n'
  exit 1
fi

# Normalize CRLF -> LF (no-op if already LF)
tr -d '\r' < "$INPUT_FILE" > "${INPUT_FILE}.unix" && mv "${INPUT_FILE}.unix" "$INPUT_FILE"

line_no=0
processed=0
skipped=0

# Use process substitution so the while loop runs in the current shell.
# tail -n +2 skips header. The "|| [ -n "$sample" ]" ensures a last line
# without newline still gets processed.
while IFS=, read -r sample forward_file reverse_file skip1 skip2 skip3 skip4 skip5 skip6 skip7 skip8 len maxdel maxins junk || [ -n "${sample:-}" ]
do
  line_no=$((line_no + 1))

  # Remove surrounding double quotes and trim leading/trailing spaces
  sample="${sample#\"}"; sample="${sample%\"}"
  sample="${sample#"${sample%%[![:space:]]*}"}"; sample="${sample%"${sample##*[![:space:]]}"}"
  forward_file="${forward_file#\"}"; forward_file="${forward_file%\"}"
  forward_file="${forward_file#"${forward_file%%[![:space:]]*}"}"; forward_file="${forward_file%"${forward_file##*[![:space:]]}"}"
  reverse_file="${reverse_file#\"}"; reverse_file="${reverse_file%\"}"
  reverse_file="${reverse_file#"${reverse_file%%[![:space:]]*}"}"; reverse_file="${reverse_file%"${reverse_file##*[![:space:]]}"}"

  # Skip completely empty lines
  if [ -z "${sample//[[:space:]]/}" ] && [ -z "${forward_file//[[:space:]]/}" ] && [ -z "${reverse_file//[[:space:]]/}" ]; then
    echo "Line $line_no: EMPTY or blank — skipping" | tee -a "$LOGFILE"
    skipped=$((skipped + 1))
    continue
  fi

  echo "Line $line_no: Processing sample '$sample' (fwd='$forward_file', rev='$reverse_file')" | tee -a "$LOGFILE"

  # Validate numeric fields (len must be integer)
  if [ -z "${len:-}" ] || [ -z "${maxdel:-}" ] || [ -z "${maxins:-}" ] || ! [[ "$len" =~ ^[0-9]+$ ]]; then
    echo "  SKIP: missing or invalid numeric fields on line $line_no: len='$len' maxdel='$maxdel' maxins='$maxins'" | tee -a "$LOGFILE"
    skipped=$((skipped + 1))
    continue
  fi

  minlen=$(( len - maxdel ))
  maxlen=$(( len + maxins ))

  forward_path="$DIR/$forward_file"
  reverse_path="$DIR/$reverse_file"
  log_path="$DIR/logs/${sample}-pandaseq-log.txt"
  merged_path="$DIR/merged/${sample}-merged.fastq"

  # Existence checks
  if [ ! -f "$forward_path" ]; then
    echo "  ERROR: forward not found: $forward_path" | tee -a "$LOGFILE"
    skipped=$((skipped + 1))
    continue
  fi
  if [ ! -f "$reverse_path" ]; then
    echo "  ERROR: reverse not found: $reverse_path" | tee -a "$LOGFILE"
    skipped=$((skipped + 1))
    continue
  fi

  # Build args array (keeps tokens separate, handles spaces)
  args=(
    -f "$forward_path"
    -r "$reverse_path"
    -g "$log_path"
    -F
    -d bFSrk
    -l "$minlen"
    -L "$maxlen"
    -w "$merged_path"
  )

  # Print & log the exact command (safely quoted)
  printf '  Running: %s ' "$PANDASEQ_BIN" | tee -a "$LOGFILE"
  for a in "${args[@]}"; do printf '%q ' "$a" | tee -a "$LOGFILE"; done
  printf '\n' | tee -a "$LOGFILE"

  # Run and capture combined stdout+stderr
  if output="$("$PANDASEQ_BIN" "${args[@]}" 2>&1)"; then
    echo "  OK: pandaseq finished for $sample" | tee -a "$LOGFILE"
    printf '%s\n' "$output" >> "$LOGFILE"

    # -------------------------
    # gzip merged FASTQ output
    # -------------------------
    if [ -s "$merged_path" ]; then
      # Use pigz if available (faster), otherwise gzip
      if command -v pigz >/dev/null 2>&1; then
        echo "  Compressing with pigz: $(basename "$merged_path")" | tee -a "$LOGFILE"
        pigz -f "$merged_path"
      else
        echo "  Compressing with gzip: $(basename "$merged_path")" | tee -a "$LOGFILE"
        gzip -f "$merged_path"
      fi
      echo "  Created: $(basename "$merged_path").gz" | tee -a "$LOGFILE"
    else
      echo "  WARNING: merged FASTQ missing or empty — not compressing" | tee -a "$LOGFILE"
    fi

    processed=$((processed + 1))
  else
    status=$?
    echo "  FAIL: pandaseq returned $status for $sample" | tee -a "$LOGFILE"
    printf '%s\n' "$output" >> "$LOGFILE"
    skipped=$((skipped + 1))
    # continue to next sample
  fi

done < <(tail -n +2 "$INPUT_FILE")

echo "Finished: $(date)" | tee -a "$LOGFILE"
echo "Summary: processed=$processed skipped=$skipped total_lines=$line_no" | tee -a "$LOGFILE"
echo "Log saved to: $LOGFILE" | tee -a "$LOGFILE"

# Pause so terminal stays open if launched via GUI
echo
read -n1 -r -p $'\nPress any key to close this window...\n'