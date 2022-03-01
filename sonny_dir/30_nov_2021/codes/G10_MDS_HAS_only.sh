#!/bin/bash

plink --bfile NWE_autosome_pruned --genome --out NWE_autosome_pruned
plink --bfile NWE_autosome_pruned --read-genome NWE_autosome_pruned.genome --cluster --mds-plot 10 --out NWE_pruned_mds
