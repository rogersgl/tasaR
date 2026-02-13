#!/usr/bin/env zsh
# Place in same folder as Input.csv and unpaired reads
# DO NOT CHANGE COLUMN ORDER in Input.csv (is crude and checks positions, not column names)

set -uo pipefail

DIR="${0:A:h}"
INPUT_FILE="${DIR}/Input.csv"
LOGFILE="${DIR}/run_pandaseq.log"

: > "$LOGFILE"
echo "Start: $(date)" >> "$LOGFILE"

mkdir -p "$DIR/logs" "$DIR/merged"

PANDASEQ_BIN="pandaseq"   # change to /usr/local/bin/pandaseq if needed
echo "Using pandaseq: $(command -v "$PANDASEQ_BIN" 2>/dev/null || echo 'NOT FOUND')" >> "$LOGFILE"

# Normalize line endings (safe no-op if already LF)
if [[ -f "$INPUT_FILE" ]]; then
  tr -d '\r' < "$INPUT_FILE" > "${INPUT_FILE}.unix" && mv "${INPUT_FILE}.unix" "$INPUT_FILE"
else
  echo "ERROR: Input file not found: $INPUT_FILE" | tee -a "$LOGFILE"
  read -k '?Press any key to close...' >/dev/tty
  exit 1
fi

line_no=0
processed=0
skipped=0

# Use process substitution so the while loop runs in current shell and variables persist.
# tail -n +2 skips header. The read ... || [[ -n $sample ]] ensures the last line
# (which might not end with a newline) still gets processed.
while IFS=, read -r sample forward_file reverse_file skip1 skip2 skip3 skip4 skip5 skip6 skip7 skip8 len maxdel maxins junk || [[ -n "${sample:-}" ]]
do
  (( line_no++ ))
  # Trim surrounding double quotes and whitespace
  sample="${sample##\"}"; sample="${sample%%\"}"
  sample="${sample## }"; sample="${sample%% }"
  forward_file="${forward_file##\"}"; forward_file="${forward_file%%\"}"
  forward_file="${forward_file## }"; forward_file="${forward_file%% }"
  reverse_file="${reverse_file##\"}"; reverse_file="${reverse_file%%\"}"
  reverse_file="${reverse_file## }"; reverse_file="${reverse_file%% }"

  # If line is empty or sample name missing, log and skip
  if [[ -z "${sample//[[:space:]]/}" && -z "${forward_file//[[:space:]]/}" && -z "${reverse_file//[[:space:]]/}" ]]; then
    echo "Line $line_no: EMPTY or blank — skipping" | tee -a "$LOGFILE"
    (( skipped++ ))
    continue
  fi

  echo "Line $line_no: Processing sample '$sample' (fwd='$forward_file', rev='$reverse_file')" | tee -a "$LOGFILE"

  # Validate numeric fields
  if [[ -z "$len" || -z "$maxdel" || -z "$maxins" || ! "$len" =~ '^[0-9]+$' ]]; then
    echo "  SKIP: missing or invalid numeric fields on line $line_no: len='$len' maxdel='$maxdel' maxins='$maxins'" | tee -a "$LOGFILE"
    (( skipped++ ))
    continue
  fi

  minlen=$(( len - maxdel ))
  maxlen=$(( len + maxins ))

  forward_path="$DIR/$forward_file"
  reverse_path="$DIR/$reverse_file"
  log_path="$DIR/logs/${sample}-pandaseq-log.txt"
  merged_path="$DIR/merged/${sample}-merged.fastq"

  # Existence checks
  if [[ ! -f "$forward_path" ]]; then
    echo "  ERROR: forward not found: $forward_path" | tee -a "$LOGFILE"
    (( skipped++ ))
    continue
  fi
  if [[ ! -f "$reverse_path" ]]; then
    echo "  ERROR: reverse not found: $reverse_path" | tee -a "$LOGFILE"
    (( skipped++ ))
    continue
  fi

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

  # Run and capture combined stdout/stderr
  if output=$("$PANDASEQ_BIN" "${args[@]}" 2>&1); then
    echo "  OK: pandaseq finished for $sample" | tee -a "$LOGFILE"
    echo "$output" >> "$LOGFILE"

    # -------------------------
    # gzip merged FASTQ output
    # -------------------------

    if [[ -s "$merged_path" ]]; then
      # Use pigz if available (faster), otherwise gzip
      if (( $+commands[pigz] )); then
        print -r -- "  Compressing with pigz: ${merged_path:t}" | tee -a -- "$LOGFILE"
        pigz -f -- "$merged_path"
      else
        print -r -- "  Compressing with gzip: ${merged_path:t}" | tee -a -- "$LOGFILE"
        gzip -f -- "$merged_path"
      fi
      print -r -- "  Created: ${merged_path:t}.gz" | tee -a -- "$LOGFILE"
    else
      print -r -- "  WARNING: merged FASTQ missing or empty — not compressing" | tee -a -- "$LOGFILE"
    fi
    (( processed++ ))
  else
    status=$?
    echo "  FAIL: pandaseq returned $status for $sample" | tee -a "$LOGFILE"
    echo "$output" >> "$LOGFILE"
    (( skipped++ ))
    # continue to next sample
  fi

done < <(tail -n +2 "$INPUT_FILE")

echo "Finished: $(date)" | tee -a "$LOGFILE"
echo "Summary: processed=$processed skipped=$skipped total_lines=$line_no" | tee -a "$LOGFILE"
echo "Log saved to: $LOGFILE" | tee -a "$LOGFILE"

# Pause so Terminal stays open after double-click
echo
read -k '?Press any key to close this window...' >/dev/tty