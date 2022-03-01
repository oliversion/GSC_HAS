#!/bin/bash
plink --bfile 1kG_MDS4 --genome --out 1kG_MDS4
plink --bfile 1kG_MDS4 --read-genome 1kG_MDS4.genome --cluster --mds-plot 10 --out 1kG_MDS4_mds
