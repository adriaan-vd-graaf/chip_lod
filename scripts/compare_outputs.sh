#!/usr/bin/env bash
# Compare every file in results/python and results/r byte-for-byte via sha256.
# Writes results/hash_check.json (file -> python hash, r hash, match; plus all_match).
# Exits non-zero if any file differs or is missing on one side. Never loosen this comparison.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="$ROOT/results/python"
RR="$ROOT/results/r"
OUT="$ROOT/results/hash_check.json"

if command -v sha256sum >/dev/null 2>&1; then
  sha() { sha256sum "$1" | cut -d' ' -f1; }
else
  sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
fi

files=$( (ls -1 "$PY" 2>/dev/null; ls -1 "$RR" 2>/dev/null) | LC_ALL=C sort -u)
if [ -z "$files" ]; then
  echo "no result files found" >&2
  exit 1
fi

all_match=true
n=0
total=$(printf '%s\n' "$files" | wc -l | tr -d ' ')
{
  printf '{\n  "files": {\n'
  for f in $files; do
    n=$((n + 1))
    hp=null; hr=null
    [ -f "$PY/$f" ] && hp="\"$(sha "$PY/$f")\""
    [ -f "$RR/$f" ] && hr="\"$(sha "$RR/$f")\""
    if [ "$hp" != null ] && [ "$hp" = "$hr" ]; then match=true; else match=false; all_match=false; fi
    sep=","; [ "$n" -eq "$total" ] && sep=""
    printf '    "%s": {\n      "python": %s,\n      "r": %s,\n      "match": %s\n    }%s\n' "$f" "$hp" "$hr" "$match" "$sep"
    printf '%-32s %s\n' "$f" "$([ "$match" = true ] && echo MATCH || echo DIFFER)" >&2
  done
  printf '  },\n  "n_files": %d,\n  "all_match": %s\n}\n' "$total" "$all_match"
} > "$OUT"

echo "all_match: $all_match (written to results/hash_check.json)" >&2
[ "$all_match" = true ]
