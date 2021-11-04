#!/usr/bin/env bash
# -*- coding: utf-8 -*-
# 
# Perform transformations of label files and extracts 
#   the mean timeseries from each label/atlas ROI.


#######################################
# Prints usage to the command line interface.
# Arguments:
#   None
#######################################
Usage(){
  cat << USAGE

  Usage: 
      
      $(basename ${0}) -f <nifti> -o <nifti> --label <nifti> --lin-xfm <mat> [--no-cleanup]
      $(basename ${0}) -f <nifti> -o <nifti> --label <nifti> --nl-xfm <warp> [--no-cleanup]

  Required arguments
    -f, --func                      Input 4D functional image file
    -o, --out-prefix                Output prefix name
    --label                         Input atlas label file

  Mutually exclusive arguments
    --lin-xfm                       Input LINEAR transformation matrix file for the standard -> native (func) transform
    --nl-xfm                        Input NON-LINEAR transformation warp file for the standard -> native (func) transform
  
  Optional arguments
    --no-cleanup                    DO NOT perform clean-up of temporary working directories
    -h, -help, --help               Prints the help menu, then exits.

USAGE
  exit 1
}


#######################################
# Prints message to the command line interface
#   in some arbitrary color.
# Arguments:
#   msg
#######################################
echo_color(){
  msg='\033[0;'"${@}"'\033[0m'
  echo -e ${msg} 
}


#######################################
# Prints message to the command line interface
#   in red.
# Arguments:
#   msg
#######################################
echo_red(){
  echo_color '31m'"${@}"
}


#######################################
# Prints message to the command line interface
#   in green.
# Arguments:
#   msg
#######################################
echo_green(){
  echo_color '32m'"${@}"
}


#######################################
# Prints message to the command line interface
#   in blue.
# Arguments:
#   msg
#######################################
echo_blue(){
  echo_color '36m'"${@}"
}


#######################################
# Prints message to the command line interface
#   in red when an error is intened to be raised.
# Arguments:
#   msg
#######################################
exit_error(){
  echo_red "${@}"
  exit 1
}


#######################################
# Logs the command to file, and executes (runs) the command.
# Globals:
#   log
#   err
# Arguments:
#   Command to be logged and performed.
#######################################
run(){
  echo "${@}"
  "${@}" >>${log} 2>>${err}
  if [[ ! ${?} -eq 0 ]]; then
    echo "failed: see log files ${log} ${err} for details"
    exit 1
  fi
  echo "-----------------------"
}


#######################################
# Performs LINEAR transformations of input atlas label files,
#   transforming the labels from atlas (standard) space to 
#   functional space.
# Arguments:
#   -f, --func: Input functional 4D image file
#   -o, --out-prefix: Output name prefix
#   --lin-xfm: Linear transformation matrix (.mat) file
#   --label: Atlas (standard) labels
#   --no-cleanup: Do NOT perform clean-up of the temporary working directory
# Outputs:
#   label file in native functional space
#######################################
lin_xfm(){
  # Set defaults
  cleanup="true"

  # Parse arguments
  while [[ ${#} -gt 0 ]]; do
    case "${1}" in
      -f|--func) shift; local func=${1} ;;
      -o|--out-prefix) shift; local out_prefix=${1} ;;
      --lin-xfm) shift; local lin_xfm=${1} ;;
      --label) shift; local label=${1} ;;
      --no-cleanup) local cleanup="false" ;;
      -*) echo_red "$(basename ${0}): Unrecognized option ${1}" >&2; Usage; ;;
      *) break ;;
    esac
    shift
  done

  local outdir=$(dirname ${out_prefix})
  local tmpdir=${outdir}/tmp_${RANDOM}

  # Create working directory
  if [[ ! -d ${tmpdir} ]]; then
    run mkdir -p ${tmpdir}
  fi

  # Set current directory, change to working directory
  local cwd=$(pwd)
  cd ${tmpdir}

  # Create reference image
  run fslroi ${func} func0 0 1

  # Apply linear xfm to label file
  run flirt -init ${lin_xfm} -in ${label} -ref func0 -applyxfm -out ${out_prefix} -interp nearestneighbour

  # Clean-up
  cd ${cwd}

  if [[ "${cleanup}" == "true" ]]; then
    run rm -rf ${tmpdir}
  fi
}


#######################################
# Performs NON-LINEAR transformations of input atlas label files,
#   transforming the labels from atlas (standard) space to 
#   functional space.
# Arguments:
#   -f, --func: Input functional 4D image file
#   -o, --out-prefix: Output name prefix
#   --nl-xfm: Non-linear transformation warp image file
#   --label: Atlas (standard) labels
#   --no-cleanup: Do NOT perform clean-up of the temporary working directory
# Outputs:
#   label file in native functional space
#######################################
nl_xfm(){
  # Set defaults
  cleanup="true"

  # Parse arguments
  while [[ ${#} -gt 0 ]]; do
    case "${1}" in
      -f|--func) shift; local func=${1} ;;
      -o|--out-prefix) shift; local out_prefix=${1} ;;
      --nl-xfm) shift; local nl_xfm=${1} ;;
      --label) shift; local label=${1} ;;
      --no-cleanup) local cleanup="false" ;;
      -*) echo_red "$(basename ${0}): Unrecognized option ${1}" >&2; Usage; ;;
      *) break ;;
    esac
    shift
  done

  local outdir=$(dirname ${out_prefix})
  local tmpdir=${outdir}/tmp_${RANDOM}

  # Create working directory
  if [[ ! -d ${tmpdir} ]]; then
    run mkdir -p ${tmpdir}
  fi

  # Set current directory, change to working directory
  local cwd=$(pwd)
  cd ${tmpdir}

  # Create reference image
  run fslroi ${func} func0 0 1

  # Apply non-linear xfm to label file
  run applywarp --in=${label} --ref=func0 --out=${out_prefix} --interp=nn --warp=${nl_xfm} --abs

  # Clean-up
  cd ${cwd}

  if [[ "${cleanup}" == "true" ]]; then
    run rm -rf ${tmpdir}
  fi
}


#######################################
# Main function that parses arguments and executes a
#   series of commands to:
#     * Transform the atlas label file to functional native space
#     * Extract the mean timeseries of each ROI
#     * Compute the adjacency/correlation matrix.
# Arguments:
#   -f, --func: Input functional 4D image file
#   -o, --out-prefix: Output name prefix
#   --lin-xfm: Linear transformation matrix (.mat) file
#   --nl-xfm: Non-linear transformation warp image file
#   --label: Atlas (standard) labels
#   --no-cleanup: Do NOT perform clean-up of the temporary working directory
# Outputs:
#   label file in native functional space
#######################################
main(){
  #
  # Parse arguments
  #============================

  # Check arguments
  if [[ ${#} -lt 1 ]]; then
    Usage >&2
    exit 1
  fi

  # Set defaults
  local cleanup="true"

  while [[ ${#} -gt 0 ]]; do
    case "${1}" in
      -f|--func) shift; local func=${1} ;;
      -o|--out-prefix) shift; local out_prefix=${1} ;;
      --lin-xfm) shift; local lin_xfm=${1} ;;
      --nl-xfm) shift; local nl_xfm=${1} ;;
      --label) shift; local label=${1} ;;
      --no-cleanup) local cleanup="false" ;;
      -h|-help|--help) Usage; ;;
      -*) echo_red "$(basename ${0}): Unrecognized option ${1}" >&2; Usage; ;;
      *) break ;;
    esac
    shift
  done

  #
  # Dependency checks
  #============================

  if ! hash flirt 2>/dev/null; then
    exit_error "FSL is not installed or added to the system path. Please check. Exiting..."
  fi

  #
  # Check arguments
  #============================

  if [[ -z ${func} ]] || [[ ! -f ${func} ]]; then
    exit_error "Functional file does not exist or was not specified."
  fi

  if [[ -z ${label} ]] || [[ ! -f ${label} ]]; then
    exit_error "Atlas label file does not exist or was not specified."
  fi

  if [[ -z ${out_prefix} ]]; then
    exit_error "Output prefix was not specified."
  fi

  if [[ -z ${lin_xfm} ]] && [[ -z ${nl_xfm} ]]; then
    exit_error "Neither linear nor non-linear transformations were specified."
  elif [[ ! -z ${lin_xfm} ]] && [[ ! -z ${nl_xfm} ]]; then
    exit_error "Both linear and non-linear transformations were specified. These options are mutually exclusive."
  elif [[ ! -z ${lin_xfm} ]] && [[ ! -f ${lin_xfm} ]]; then
    exit_error "Input linear transfromation matrix file does not exist."
  elif [[ ! -z ${nl_xfm} ]] && [[ ! -f ${nl_xfm} ]]; then
    exit_error "Input non-linear transfromation warp file does not exist."
  fi

  #
  # Define output directory
  #============================

  # Define log files
  local outdir=$(dirname ${out_prefix})

  # Create output directory
  if [[ ! -d ${outdir} ]]; then
    mkdir -p ${outdir}
  fi

  local outdir=$(realpath ${outdir})

  log=${outdir}/gen_ts_mat.log
  err=${outdir}/gen_ts_mat.err

  #
  # Perform transforms
  #============================

  if [[ "${cleanup}" == "true" ]]; then
    local clean=""
  elif [[ "${cleanup}" == "false" ]]; then
    local clean="--no-cleanup"
  fi

  if [[ -f ${lin_xfm} ]]; then
    run lin_xfm --func ${func} --out-prefix ${out_prefix} --lin-xfm ${lin_xfm} --label ${label} ${clean}
  elif [[ -f ${nl_xfm} ]]; then
    run nl_xfm --func ${func} --out-prefix ${out_prefix} --nl-xfm ${nl_xfm} --label ${label} ${clean}
  fi

  #
  # Extract mean timeseries
  #============================

  run fslmeants -i ${func} -o ${out_prefix}_mean_ts_rois_mat.txt --label=${out_prefix}.nii.gz

  #
  # Compute correlation matrix
  #============================

  prog_dir=$(dirname $(realpath ${0}))
  run ${prog_dir}/corr_mat.py --matrix-file=${out_prefix}_mean_ts_rois_mat.txt --output=${out_prefix}_corr_mat.txt
  run ${prog_dir}/corr_mat.py --matrix-file=${out_prefix}_mean_ts_rois_mat.txt --output=${out_prefix}_corr_mat_Z_xfm.txt -z

  #
  # Clean-up
  #============================
  
  if [[ "${cleanup}" == "true" ]]; then
    run rm ${out_prefix}.nii.gz # ${out_prefix}_mean_rois_mat.txt
  fi
  
  # Successful exit status
  exit 0
}

# Main function
main "${@}"
