#! /bin/bash

chemsh_PATH=/home/uabqut17/soft/chemshell-3.7/scripts/chemsh

PDB=$1
C=`echo "$2" | cut -d'.' -f1`


echo "read_pdb file=$PDB coords=dummy.coords" > chemsh.tmp
echo "write_pdb file=${C}.pdb coords=${C}.c" >> chemsh.tmp

$chemsh_PATH chemsh.tmp

rm chemsh.tmp
