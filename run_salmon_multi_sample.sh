#!/usr/bin/env bash
set -euo pipefail

# Multi-sample Salmon quant pipeline for macOS + Jupyter notebook usage.
# Samples: PC1, PC2, PC3, PC7, PC8, PC9
#
# Output naming convention:
#   quant_PC1_regular.sf
#   quant_PC1_bootstrap50.sf
#   meta_PC1_regular.json
#   salmon_PC1_regular.log
#
# Run as:
#   bash "run_salmon_multi_sample.sh"
# or paste into a %%bash notebook cell.

THREADS="${THREADS:-8}"
SALMON_VERSION="${SALMON_VERSION:-1.11.4}"
RUN_MODE="${RUN_MODE:-regular}"          # regular | bootstrap50 | both
BOOTSTRAPS="${BOOTSTRAPS:-0}"            # overrides bootstrap count when >0
GIBBS_SAMPLES="${GIBBS_SAMPLES:-0}"      # optional future extension; guarded against bootstraps
WRITE_UNMAPPED_NAMES="${WRITE_UNMAPPED_NAMES:-0}"

SAMPLES=(PC1 PC2 PC3 PC7 PC8 PC9)

# Resolve project directory whether script is launched from:
# 1) project root ("Pear Irradiation"), or
# 2) parent directory containing "Pear Irradiation".
if [[ -d "Pear Irradiation" ]]; then
  PROJECT_DIR="$(cd "Pear Irradiation" && pwd)"
else
  PROJECT_DIR="$(pwd)"
fi

TRANSCRIPTOME="${PROJECT_DIR}/Bartlett_v2.0.cds.fasta"
READS_DIR="${PROJECT_DIR}/Paired_FQ_Files"
INDEX_DIR="${PROJECT_DIR}/bartlett_index"
QUANTS_DIR="${PROJECT_DIR}/quants"
LOG_DIR="${PROJECT_DIR}/logs"
INSTALL_DIR="${PROJECT_DIR}/tools"
LOCAL_BIN="${INSTALL_DIR}/bin"

MAIN_LOG="${LOG_DIR}/pipeline_$(date +%Y%m%d_%H%M%S).log"
COMPARISON_TSV="${LOG_DIR}/sample_run_comparison.tsv"
UNMAPPED_TSV="${LOG_DIR}/unmapped_summary.tsv"

mkdir -p "${LOG_DIR}" "${QUANTS_DIR}" "${LOCAL_BIN}"
touch "${MAIN_LOG}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "${MAIN_LOG}"
}

die() {
  log "ERROR: $*"
  exit 1
}

download_salmon_binary() {
  local arch tarball url candidate d
  arch="$(uname -m)"

  case "${arch}" in
    arm64) tarball="salmon-${SALMON_VERSION}_macos_arm64.tar.gz" ;;
    x86_64) tarball="salmon-${SALMON_VERSION}_macos_x86_64.tar.gz" ;;
    *) die "Unsupported architecture '${arch}' for direct Salmon binary install." ;;
  esac

  url="https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/${tarball}"
  log "Downloading Salmon binary from ${url}"
  curl -fL --retry 3 --retry-delay 3 -o "${INSTALL_DIR}/${tarball}" "${url}" || return 1

  log "Extracting ${tarball}"
  tar -xzf "${INSTALL_DIR}/${tarball}" -C "${INSTALL_DIR}"

  candidate=""
  for d in "${INSTALL_DIR}"/*; do
    if [[ -x "${d}/bin/salmon" ]]; then
      candidate="${d}/bin/salmon"
      break
    fi
  done

  [[ -n "${candidate}" ]] || return 1
  cp "${candidate}" "${LOCAL_BIN}/salmon"
  chmod +x "${LOCAL_BIN}/salmon"
  export PATH="${LOCAL_BIN}:${PATH}"
  log "Installed Salmon binary to ${LOCAL_BIN}/salmon"
}

ensure_salmon() {
  if command -v salmon >/dev/null 2>&1; then
    log "Found Salmon at $(command -v salmon)"
    salmon --version
    return 0
  fi

  log "Salmon not found; trying Homebrew."
  if command -v brew >/dev/null 2>&1; then
    if brew install salmon && command -v salmon >/dev/null 2>&1; then
      salmon --version
      return 0
    fi
    log "Homebrew method failed; trying next method."
  else
    log "Homebrew not available."
  fi

  # Local bundled Salmon binaries if present.
  for candidate in \
    "${PROJECT_DIR}/salmon-macos-arm64/bin/salmon" \
    "${PROJECT_DIR}/salmon-linux-x86_64/bin/salmon"; do
    if [[ -x "${candidate}" ]]; then
      cp "${candidate}" "${LOCAL_BIN}/salmon"
      chmod +x "${LOCAL_BIN}/salmon"
      export PATH="${LOCAL_BIN}:${PATH}"
      if salmon --version; then
        log "Using bundled Salmon binary: ${candidate}"
        return 0
      fi
    fi
  done

  log "Trying direct binary install from Salmon releases."
  if download_salmon_binary && salmon --version; then
    return 0
  fi

  # Docker fallback if available.
  if command -v docker >/dev/null 2>&1; then
    log "Trying Docker fallback."
    docker pull combinelab/salmon || true
    if docker run --rm combinelab/salmon salmon --version; then
      export SALMON_DOCKER_MODE="1"
      return 0
    fi
  fi

  die "Unable to obtain a working Salmon installation."
}

salmon_cmd() {
  if [[ "${SALMON_DOCKER_MODE:-0}" == "1" ]]; then
    docker run --rm \
      -v "${PROJECT_DIR}:${PROJECT_DIR}" \
      -w "${PROJECT_DIR}" \
      combinelab/salmon salmon "$@"
  else
    salmon "$@"
  fi
}

validate_inputs() {
  [[ -f "${TRANSCRIPTOME}" ]] || die "Transcriptome FASTA missing: ${TRANSCRIPTOME}"
  [[ -d "${READS_DIR}" ]] || die "Reads directory missing: ${READS_DIR}"
}

find_sample_pair() {
  local sample="$1"
  local r1_candidates=(
    "${READS_DIR}/${sample}_R1"*.fq
    "${READS_DIR}/${sample}_R1"*.fastq
    "${READS_DIR}/${sample}_R1"*.fq.gz
    "${READS_DIR}/${sample}_R1"*.fastq.gz
  )
  local r1=""
  local r2=""
  local p

  shopt -s nullglob
  for p in "${r1_candidates[@]}"; do
    if [[ -f "${p}" ]]; then
      r1="${p}"
      break
    fi
  done
  shopt -u nullglob

  [[ -n "${r1}" ]] || die "No ${sample} R1 file found in ${READS_DIR} (.fq/.fastq/.gz supported)."

  r2="${r1/_R1/_R2}"
  [[ -f "${r2}" ]] || die "Missing matching R2 for ${r1}. Expected ${r2}"

  SAMPLE_R1="${r1}"
  SAMPLE_R2="${r2}"
}

build_index_if_needed() {
  if [[ -d "${INDEX_DIR}" && -f "${INDEX_DIR}/info.json" && -f "${INDEX_DIR}/versionInfo.json" ]]; then
    log "Index appears valid. Skipping rebuild: ${INDEX_DIR}"
    return 0
  fi

  log "Building Salmon index: ${INDEX_DIR}"
  salmon_cmd index -t "${TRANSCRIPTOME}" -i "${INDEX_DIR}" -k 31
  [[ -f "${INDEX_DIR}/info.json" ]] || die "Index build failed: info.json not found."
}

bootstraps_for_mode() {
  local mode="$1"
  case "${mode}" in
    regular) echo "0" ;;
    bootstrap50) echo "50" ;;
    *) die "Unsupported internal mode '${mode}'." ;;
  esac
}

trial_label_for_bootstraps() {
  local n="$1"
  if [[ "${n}" -eq 0 ]]; then
    echo "regular"
  else
    echo "bootstrap${n}"
  fi
}

extract_meta_field() {
  local meta_json="$1"
  local field="$2"
  python3 - "${meta_json}" "${field}" <<'PY'
import json, sys
meta_path, field = sys.argv[1], sys.argv[2]
with open(meta_path, "r", encoding="utf-8") as fh:
    data = json.load(fh)
print(data.get(field, ""))
PY
}

append_metrics_row() {
  local sample="$1"
  local trial="$2"
  local out_dir="$3"
  local meta_json="${out_dir}/aux_info/meta_info.json"
  local num_processed num_mapped percent_mapped num_bootstraps num_unmapped percent_unmapped
  local row

  [[ -f "${meta_json}" ]] || die "Missing Salmon metadata file: ${meta_json}"

  num_processed="$(extract_meta_field "${meta_json}" "num_processed")"
  num_mapped="$(extract_meta_field "${meta_json}" "num_mapped")"
  percent_mapped="$(extract_meta_field "${meta_json}" "percent_mapped")"
  num_bootstraps="$(extract_meta_field "${meta_json}" "num_bootstraps")"

  [[ -n "${num_processed}" && -n "${num_mapped}" && -n "${percent_mapped}" && -n "${num_bootstraps}" ]] \
    || die "Required mapping fields missing in ${meta_json}"

  num_unmapped=$((num_processed - num_mapped))
  percent_unmapped="$(python3 - "${percent_mapped}" <<'PY'
import sys
v = float(sys.argv[1])
print(f"{100.0 - v:.6f}")
PY
)"

  row=$(printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' \
    "${sample}" "${trial}" "${num_processed}" "${num_mapped}" "${num_unmapped}" \
    "${percent_mapped}" "${percent_unmapped}" "${num_bootstraps}")

  printf '%s\n' "${row}" >> "${COMPARISON_TSV}"
  printf '%s\n' "${row}" >> "${UNMAPPED_TSV}"

  log "Metrics ${sample}/${trial}: processed=${num_processed}, mapped=${num_mapped}, unmapped=${num_unmapped}, percent_mapped=${percent_mapped}, percent_unmapped=${percent_unmapped}, num_bootstraps=${num_bootstraps}"
}

tag_run_outputs() {
  local sample="$1"
  local trial="$2"
  local out_dir="$3"
  local tag="${sample}_${trial}"
  local file rel base dir new_base new_path

  find "${out_dir}" -type f -print0 | while IFS= read -r -d '' file; do
    rel="${file#${out_dir}/}"
    base="$(basename "${rel}")"
    dir="$(dirname "${rel}")"

    case "${base}" in
      *_"${tag}".*|*_"${tag}") continue ;;
    esac

    if [[ "${base}" == *.* ]]; then
      new_base="${base%.*}_${tag}.${base##*.}"
    else
      new_base="${base}_${tag}"
    fi

    if [[ "${dir}" == "." ]]; then
      new_path="${out_dir}/${new_base}"
    else
      new_path="${out_dir}/${dir}/${new_base}"
    fi

    mv "${file}" "${new_path}"
  done
}

run_quant_sample_trial() {
  local sample="$1"
  local bootstraps="$2"
  local r1="$3"
  local r2="$4"
  local trial
  local out_dir
  local salmon_log
  local control_log
  local -a quant_extra=()
  local -a quant_cmd=()
  local unmapped_count

  trial="$(trial_label_for_bootstraps "${bootstraps}")"
  out_dir="${QUANTS_DIR}/${sample}_${trial}"
  salmon_log="${LOG_DIR}/salmon_${sample}_${trial}.log"
  control_log="${LOG_DIR}/run_${sample}_${trial}.log"

  rm -rf "${out_dir}"
  mkdir -p "${out_dir}"

  if [[ "${bootstraps}" -gt 0 ]]; then
    quant_extra+=(--numBootstraps "${bootstraps}")
  fi
  if [[ "${GIBBS_SAMPLES}" -gt 0 ]]; then
    quant_extra+=(--numGibbsSamples "${GIBBS_SAMPLES}")
  fi
  if [[ "${WRITE_UNMAPPED_NAMES}" == "1" ]]; then
    quant_extra+=(--writeUnmappedNames)
  fi

  log "Processing ${sample} trial=${trial} out=${out_dir}"
  printf 'Processing sample: %s (%s)\n' "${sample}" "${trial}" | tee -a "${LOG_DIR}/quant_run.log" "${control_log}"

  quant_cmd=(
    quant
    -i "${INDEX_DIR}"
    -l A
    -1 "${r1}"
    -2 "${r2}"
    -p "${THREADS}"
    --validateMappings
  )

  if [[ "${#quant_extra[@]}" -gt 0 ]]; then
    quant_cmd+=("${quant_extra[@]}")
  fi
  quant_cmd+=(-o "${out_dir}")

  if salmon_cmd "${quant_cmd[@]}" > "${salmon_log}" 2>&1; then
    printf 'SUCCESS %s (%s)\n' "${sample}" "${trial}" | tee -a "${LOG_DIR}/quant_run.log" "${control_log}"
  else
    printf 'FAILED %s (%s) (see %s)\n' "${sample}" "${trial}" "${salmon_log}" | tee -a "${LOG_DIR}/quant_failures.log" "${control_log}"
    die "Quantification failed for ${sample} trial=${trial}."
  fi

  [[ -s "${out_dir}/quant.sf" ]] || die "Missing or empty ${out_dir}/quant.sf"

  if [[ "${WRITE_UNMAPPED_NAMES}" == "1" ]]; then
    [[ -f "${out_dir}/aux_info/unmapped_names.txt" ]] || die "Expected unmapped names file missing for ${sample} trial=${trial}"
    unmapped_count="$(python3 - "${out_dir}/aux_info/unmapped_names.txt" <<'PY'
import sys
count = 0
with open(sys.argv[1], "r", encoding="utf-8") as fh:
    for _ in fh:
        count += 1
print(count)
PY
)"
    log "${sample} trial=${trial} unmapped_names_count=${unmapped_count}"
  fi

  append_metrics_row "${sample}" "${trial}" "${out_dir}"
  tag_run_outputs "${sample}" "${trial}" "${out_dir}"

  [[ -s "${out_dir}/quant_${sample}_${trial}.sf" ]] || die "Missing tagged file ${out_dir}/quant_${sample}_${trial}.sf"
  [[ -f "${out_dir}/aux_info/meta_info_${sample}_${trial}.json" ]] || die "Missing tagged meta file for ${sample} trial=${trial}"

  if [[ "${WRITE_UNMAPPED_NAMES}" == "1" ]]; then
    [[ -f "${out_dir}/aux_info/unmapped_names_${sample}_${trial}.txt" ]] || die "Missing tagged unmapped names file for ${sample} trial=${trial}"
  fi

  log "Tagged outputs created for ${sample}/${trial}"
}

prepare_summary_files() {
  : > "${LOG_DIR}/quant_run.log"
  : > "${LOG_DIR}/quant_failures.log"
  salmon_cmd --version | tee "${LOG_DIR}/salmon_version.txt"

  printf 'sample\ttrial\tnum_processed\tnum_mapped\tnum_unmapped\tpercent_mapped\tpercent_unmapped\tnum_bootstraps\n' > "${COMPARISON_TSV}"
  printf 'sample\ttrial\tnum_processed\tnum_mapped\tnum_unmapped\tpercent_mapped\tpercent_unmapped\tnum_bootstraps\n' > "${UNMAPPED_TSV}"
}

run_all_samples_and_trials() {
  local modes=()
  local sample mode bootstraps
  local processed=0
  local failed=0

  case "${RUN_MODE}" in
    regular|baseline) modes=("regular") ;;
    bootstrap50) modes=("bootstrap50") ;;
    both) modes=("regular" "bootstrap50") ;;
    *) die "RUN_MODE must be regular, bootstrap50, or both. Got '${RUN_MODE}'." ;;
  esac

  for sample in "${SAMPLES[@]}"; do
    find_sample_pair "${sample}"
    log "Detected input pair for ${sample}: ${SAMPLE_R1} | ${SAMPLE_R2}"
  done

  for sample in "${SAMPLES[@]}"; do
    for mode in "${modes[@]}"; do
      bootstraps="$(bootstraps_for_mode "${mode}")"
      if [[ "${mode}" == "bootstrap50" && "${BOOTSTRAPS}" -gt 0 ]]; then
        bootstraps="${BOOTSTRAPS}"
      fi

      run_quant_sample_trial "${sample}" "${bootstraps}" "${SAMPLE_R1}" "${SAMPLE_R2}"
      processed=$((processed + 1))
    done
  done

  printf 'processed=%s\nfailed=%s\n' "${processed}" "${failed}" | tee "${LOG_DIR}/summary.txt"
}

validate_requested_sampling() {
  [[ "${BOOTSTRAPS}" =~ ^[0-9]+$ ]] || die "BOOTSTRAPS must be a non-negative integer."
  [[ "${GIBBS_SAMPLES}" =~ ^[0-9]+$ ]] || die "GIBBS_SAMPLES must be a non-negative integer."
  if [[ "${BOOTSTRAPS}" -gt 0 && "${GIBBS_SAMPLES}" -gt 0 ]]; then
    die "BOOTSTRAPS and GIBBS_SAMPLES cannot both be > 0 in a single run."
  fi
}

main() {
  log "Starting multi-sample Salmon pipeline (samples=${SAMPLES[*]}, RUN_MODE=${RUN_MODE}, BOOTSTRAPS=${BOOTSTRAPS}, WRITE_UNMAPPED_NAMES=${WRITE_UNMAPPED_NAMES})."
  validate_inputs
  validate_requested_sampling
  ensure_salmon
  build_index_if_needed
  prepare_summary_files
  run_all_samples_and_trials
  log "Pipeline completed successfully. Logs: ${LOG_DIR} | Outputs under: ${QUANTS_DIR}"
}

main "$@"