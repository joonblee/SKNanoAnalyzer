#!/bin/bash
# Source-safe entry point: all settings and exits stay in this subshell.
(

##############
### setups ###
##############

# run setups
RUN_DT=true
RUN_MC=true
RUN_QCDonly=false   # Available only with 'RUN_MC=true'
RUN_SIG=true

UseSkim=true
# NOTE: TriggerEfficiency requires TrigObj branches in the input.
#       If UseSkim=true, the Skim_NIsoMuon output must retain TrigObj_* branches.

# analysis setups
analysis="NIsoMuon"
skim="Skim_NIsoMuon"
TriggerSets=("HighPtMuon")

signalset_Run2="SampleLists/Run2signal.txt"
signalset_Run3="SampleLists/Run3signal.txt"

mcset_Run2="SampleLists/Run2mc.txt"
#mcset_Run3="SampleLists/Run3mc_all.txt"
mcset_Run3="SampleLists/Run3mc.txt"

qcdset_Run2="SampleLists/Run2qcd.txt"
qcdset_Run3="SampleLists/Run3qcd.txt"

systset_Run2="SampleLists/Run2Syst.txt"
systset_Run3="SampleLists/Run3Syst.txt"

effset_Run3="SampleLists/Run3eff.txt"


nBatch=50
batchname=""       # e.g. "NIsoMuon"
extra_args=()      # e.g. ("--no_exec") or ("--memory" "4000")


# CLI options are local to the subshell, including when this file is sourced.
era_selection="Run3"
flag_selection="nominal"
top_pt_reweight="${TOP_PT_REWEIGHT:-false}"
dry_run=false
print_usage() {
  cat <<'HELP'
Usage: source run.sh [options]     (bash run.sh also works)
  --era, --eras LIST   Run3 (default), ERA, or comma-separated eras/runs
  --flags LIST        nominal (default), RunSyst, RunXSecSyst, or all
                      Comma-separated modes are submitted independently.
                      MuonIDEfficiency and TriggerEfficiency are also accepted.
  --top-pt-reweight [true|false]   Optional ttbar-only study (default: false)
  --no-top-pt-reweight            Explicitly disable the study
  --dry-run                      Print commands without building/submitting
  -h, --help
Examples:
  source run.sh
  source run.sh --era 2022EE --flags nominal,RunSyst
  source run.sh --era Run3 --flags all --top-pt-reweight
  TOP_PT_REWEIGHT=true bash run.sh --flags all
HELP
}
while (( $# )); do
  case "$1" in
    --era|--eras|--flags)
      if (( $# < 2 )) || [[ "$2" == --* ]]; then echo "[ERROR] $1 requires a value" >&2; exit 2; fi
      if [[ "$1" == "--flags" ]]; then flag_selection="$2"; else era_selection="$2"; fi
      shift 2 ;;
    --top-pt-reweight)
      top_pt_reweight=true
      if [[ "${2:-}" == "true" || "${2:-}" == "false" ]]; then top_pt_reweight="$2"; shift; fi
      shift ;;
    --no-top-pt-reweight) top_pt_reweight=false; shift ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help) print_usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; print_usage; exit 2 ;;
  esac
done
case "$top_pt_reweight" in true|false) ;; *) echo "[ERROR] TOP_PT_REWEIGHT must be true or false" >&2; exit 2 ;; esac

Eras=()
IFS=',' read -r -a selections <<< "$era_selection"
for selection in "${selections[@]}"; do
  case "$selection" in
    Run2|run2) expanded=(2016preVFP 2016postVFP 2017 2018) ;;
    Run3|run3) expanded=(2022 2022EE 2023 2023BPix) ;;
    Run2+3|run2+3|all) expanded=(2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix) ;;
    2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix) expanded=("$selection") ;;
    *) echo "[ERROR] Unknown era: $selection" >&2; exit 2 ;;
  esac
  for era in "${expanded[@]}"; do
    if [[ " ${Eras[*]} " != *" $era "* ]]; then Eras+=("$era"); fi
  done
done
flags=()
IFS=',' read -r -a selections <<< "$flag_selection"
for selection in "${selections[@]}"; do
  case "$selection" in
    nominal|Nominal) expanded=("") ;;
    all) expanded=("" "RunSyst" "RunXSecSyst") ;;
    RunSyst|RunXSecSyst|MuonIDEfficiency|TriggerEfficiency) expanded=("$selection") ;;
    *) echo "[ERROR] Unknown flag: $selection" >&2; exit 2 ;;
  esac
  for flag in "${expanded[@]}"; do
    if $top_pt_reweight; then
      case "$flag" in
        ""|RunSyst|RunXSecSyst) flag="TopPtReweight${flag:+,$flag}" ;;
        *) echo "[ERROR] TopPtReweight supports nominal, RunSyst and RunXSecSyst only" >&2; exit 2 ;;
      esac
    fi
    duplicate=false
    for existing in "${flags[@]}"; do [[ "$existing" == "$flag" ]] && duplicate=true; done
    $duplicate || flags+=("$flag")
  done
done
if (( ${#Eras[@]} == 0 || ${#flags[@]} == 0 )); then echo "[ERROR] Empty era/mode selection" >&2; exit 2; fi

submit_job() {
  local logfile="$1"
  shift
  if $dry_run; then
    printf '[dry-run] '
    printf '%q ' "$@"
    printf '\n'
  else
    "$@" &> "$logfile"
  fi
}

if ! $dry_run; then
  mkdir -p log
fi

list_to_input() {
  local listfile="$1"
  local prefix="$2"

  awk -v prefix="$prefix" '
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      print prefix $0
    }
  ' "$listfile" | paste -sd, -
}

prefix_input() {
  local samples="$1"
  local prefix="$2"
  local out=""
  local sample

  IFS=',' read -ra sample_array <<< "$samples"
  for sample in "${sample_array[@]}"; do
    if [[ -n "$out" ]]; then
      out+=","
    fi
    out+="${prefix}${sample}"
  done
  echo "$out"
}

is_efficiency_mode() {
  local flag="$1"
  [[ "$flag" == "MuonIDEfficiency" || "$flag" == "TriggerEfficiency" ]]
}

is_systematic_mode() {
  local flag="$1"
  [[ "$flag" == "RunSyst" || "$flag" == "RunXSecSyst" ]]
}

echo ""
echo "// ------------------------------------------------------ //"
echo "// ----------------- Run SKNanoAnalyzer ----------------- //"
echo "// ------------------------------------------------------ //"
echo ""
echo -n "// time: "
date
echo ""

for trig in "${TriggerSets[@]}"; do

  for era in "${Eras[@]}"; do

    if [[ "$era" == "2016preVFP" || "$era" == "2016postVFP" || "$era" == "2017" || "$era" == "2018" ]]; then
      dataset="SingleMuon"
      signalset="$signalset_Run2"
      mcset="$mcset_Run2"
      qcdset="$qcdset_Run2"
      systset="$systset_Run2"
    elif [[ "$era" == "2022" || "$era" == "2022EE" ]]; then
      dataset="Muon"
      signalset="$signalset_Run3"
      mcset="$mcset_Run3"
      qcdset="$qcdset_Run3"
      systset="$systset_Run3"
    elif [[ "$era" == "2023" || "$era" == "2023BPix" ]]; then
      dataset="Muon0,Muon1"
      signalset="$signalset_Run3"
      mcset="$mcset_Run3"
      qcdset="$qcdset_Run3"
      systset="$systset_Run3"
    fi

    for flag in "${flags[@]}"; do

      cmd_common=(SKNano.py -a "$analysis" -e "$era" -n "$nBatch" --nmax 4000)

      # Submit only the POWHEG ttbar samples for the optional study.
      # Keep all ordinary production branches below unchanged.
      if [[ ",$flag," == *",TopPtReweight,"* ]]; then
        if $RUN_MC; then
          prefix=""
          if $UseSkim; then prefix="${skim}_"; fi
          toppt_input=$(awk -v prefix="$prefix" \
            '/^[[:space:]]*TT(JJ|LJ|LL)_powheg([_[:alnum:]]*)[[:space:]]*$/ {print prefix $1}' \
            "$mcset" | paste -sd, -)
          if [[ -n "$toppt_input" ]]; then
            if [[ -n "$batchname" ]]; then cmd_common+=(--batchname "$batchname"); fi
            submit_job "log/submit_${era}_ttbar_${trig}__${flag}.log" \
              "${cmd_common[@]}" --userflags "$flag" "${extra_args[@]}" -i "$toppt_input"
            echo "[SKNano.py] Run ttbar TopPtReweight study: $flag era: $era"
          fi
        fi
        continue
      fi

      echo " - Era: $era"
      echo " - Trigger: $trig"
      echo " - Dataset: $dataset"
      echo " - MCset: $mcset"

      if [[ -n "$flag" ]]; then
        cmd_common+=(--userflags "$flag")
        echo " - flag: $flag"
        if is_efficiency_mode "$flag"; then
          echo " - mode: central-only efficiency measurement"
        elif [[ "$flag" == "RunSyst" ]]; then
          echo " - mode: NIsoDimuon experimental systematics"
        elif [[ "$flag" == "RunXSecSyst" ]]; then
          echo " - mode: NIsoDimuon generator theory systematics"
        fi
      else
        echo " - flag: <none>"
        echo " - mode: nominal NIsoDimuon"
      fi
      echo ""

      mcset_this="$mcset"
      
      if $RUN_QCDonly; then
        mcset_this="$qcdset"
      fi
      
      if is_systematic_mode "$flag"; then
        mcset_this="$systset"
      fi
      
      if is_efficiency_mode "$flag"; then
        if [[ "$era" == "2022" || "$era" == "2022EE" || \
              "$era" == "2023" || "$era" == "2023BPix" ]]; then
          mcset_this="$effset_Run3"
        fi
      fi

      # MuonIDEfficiency and TriggerEfficiency are central-only modes.
      # They use the ordinary nominal MC list above, not the RunXSecSyst list.

      if [[ -n "$batchname" ]]; then
        cmd_common+=(--batchname "$batchname")
      fi

      if ((${#extra_args[@]})); then
        cmd_common+=("${extra_args[@]}")
      fi

      if $RUN_DT && ! is_systematic_mode "$flag"; then
        if $UseSkim; then
          data_input=$(prefix_input "$dataset" "${skim}_")
        else
          data_input="$dataset"
        fi
      fi

      if $RUN_MC; then
        if $UseSkim; then
          mc_input=$(list_to_input "$mcset_this" "${skim}_")
        else
          mc_input=$(list_to_input "$mcset_this" "")
        fi
      fi

      if $RUN_SIG && ! is_efficiency_mode "$flag" && [[ "$flag" != "RunXSecSyst" ]]; then
        if $UseSkim; then
          sig_input=$(list_to_input "$signalset" "${skim}_")
        else
          sig_input=$(list_to_input "$signalset" "")
        fi
      fi

      if $RUN_DT && ! is_systematic_mode "$flag"; then
        submit_job "log/submit_${era}_${dataset//,/_}_${trig}${flag:+__${flag}}.log" \
          "${cmd_common[@]}" \
          -i "$data_input"
        echo "[SKNano.py] Run analyzer for data: $dataset trigger: $trig era: $era"
      else
        echo "[SKNano.py] Do not make DATA samples"
      fi

      if $RUN_MC; then
        submit_job "log/submit_${era}_MC_${trig}${flag:+__${flag}}.log" \
          "${cmd_common[@]}" \
          -i "$mc_input"
        echo "[SKNano.py] Run analyzer for MC backgrounds: $mcset_this trigger: $trig era: $era"
      else
        echo "[SKNano.py] Do not make MC samples"
      fi

      if $RUN_SIG && ! is_efficiency_mode "$flag" && [[ "$flag" != "RunXSecSyst" ]]; then
        cmd_signal=("${cmd_common[@]}")
        cmd_signal[6]=1
        submit_job "log/submit_${era}_sig_${trig}${flag:+__${flag}}.log" \
          "${cmd_signal[@]}" \
          -i "$sig_input"
        echo "[SKNano.py] Run analyzer for signals: $signalset trigger: $trig era: $era"
      elif $RUN_SIG && is_efficiency_mode "$flag"; then
        echo "[SKNano.py] Skip signal samples for efficiency mode: $flag"
      else
        echo "[SKNano.py] Do not make signal samples"
      fi

    done
  done
done


)
