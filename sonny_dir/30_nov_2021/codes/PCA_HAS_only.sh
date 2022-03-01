#!/bin/bash

## Perform PCA using king on the HAS data to obtain the covariates for further analyses
king -b HAS_EUR1.bed --pca --cpus 72 --prefix HAS_EUR_final_pc
