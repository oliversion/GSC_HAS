#!/bin/bash

### Remove non-EUR population by FID

plink --bfile HAS_MDS3 --remove-fam infer_non_EUR.txt --make-bed --out HAS_EUR1

