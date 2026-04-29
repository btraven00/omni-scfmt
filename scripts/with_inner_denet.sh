#!/usr/bin/env bash
# Inner-container denet wrapper.
#
# Usage: with_inner_denet.sh <cmd> <args...>
#
# Extracts --output_dir from the downstream arg list so the inner denet
# JSON lands next to the other per-rule artifacts (omnibench-events.jsonl,
# etc.). Falls back to CWD if --output_dir is absent.
#
# The *_profile entrypoint of each module is a 2-line script that execs
# into this wrapper, so every module gets identical inner profiling with
# zero duplication.

set -euo pipefail

out_dir="."
args=("$@")
for ((i = 1; i < ${#args[@]}; i++)); do
    if [[ "${args[i]}" == "--output_dir" && $((i + 1)) -lt ${#args[@]} ]]; then
        out_dir="${args[i + 1]}"
        break
    fi
done
mkdir -p "$out_dir"

# Use the rule-ish name of the first arg as the denet file prefix so multiple
# modules writing into the same output_dir don't clobber each other.
cmd_basename="$(basename "${1%.*}")"
denet_out="${out_dir}/${cmd_basename}.inner_denet.jsonl"

exec denet --json --quiet --out "$denet_out" run -- "$@"
