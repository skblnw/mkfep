#!/bin/bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
  echo "Usage: $0 <target_directory> <subdirectory> <number_of_windows>"
  exit 1
fi

target_dir=$(realpath "$1")
subdir="$2"
num_windows="$3"

missing_list=()
free_list=()
complex_list=()
both_list=()

# Print table header
printf "%-30s | %-15s | %-15s\n" "Directory" "Free State" "Complex State"
printf "%-30s-+-%-15s-+-%-15s\n" "------------------------------" "---------------" "---------------"

for ii in "${target_dir}"/*/"${subdir}"/; do
  relative_ii=${ii##"${target_dir}/"}
  dir_name=${relative_ii%/"${subdir}"/}
  md_files=($(find "$ii" -name "md.xvg"))
  count=${#md_files[@]}
  free_count=0
  complex_count=0
  for file in "${md_files[@]}"; do
    if [[ $file =~ "free" ]]; then
      ((free_count++))
    elif [[ $file =~ "complex" ]]; then
      ((complex_count++))
    fi
  done

  free_status="Incomplete"
  complex_status="Incomplete"

  if [ "$count" -ne $((num_windows * 2)) ]; then
    missing_list+=("$dir_name")

    if [ "$free_count" -eq "$num_windows" ]; then
      free_list+=("$dir_name")
      free_status="Complete"
    fi
    if [ "$complex_count" -eq "$num_windows" ]; then
      complex_list+=("$dir_name")
      complex_status="Complete"
    fi
  else
    both_list+=("$dir_name")
    free_status="Complete"
    complex_status="Complete"
  fi

  # Print table row
  printf "%-30s | %-15s | %-15s\n" "$dir_name" "$free_status" "$complex_status"
done
