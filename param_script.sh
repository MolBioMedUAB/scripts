#!/bin/bash

ml Amber

while getopts "m:n:" opt; do
  case $opt in
    m)
      if [ "$OPTARG" == "1" ] || [ "$OPTARG" == "2" ]; then
        mode="$OPTARG"
      else
        echo "Invalid mode: $OPTARG. Please use 1 or 2."
        exit 1
      fi
      ;;
    n)
      name="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND - 1))

if [ -z "$mode" ] || [ -z "$name" ]; then
  echo "Usage: $0 -m [1 or 2] -n [name]"
  exit 1
fi

  case $mode in
    1)
      antechamber -dr no -i "${name}.pdb" -fi pdb -o "${name}.com" -fo gcrt -gv 1
      echo "Sent *com file to Gaussian"
      ;;
    2)
      espgen -i "${name}.gesp" -o "${name}.esp"
      antechamber -i "${name}.log" -fi gout -o "${name}.ac" -fo ac -c esp -cf "${name}.esp"
      parmchk2 -i "${name}.ac" -f ac -o "${name}.frcmod" -a y -w y -s 2
      antechamber -i "${name}.ac" -fi ac -o "${name}.prepi" -fo prepi
      antechamber -i "${name}.ac" -fi ac -o "${name}.mol2" -fo mol2

      ;;
  esac
