#!/bin/bash
#$ -M mbonas@nd.edu
#$ -m ae
#$ -q long
#$ -pe smp 1
#$ -l h=(d12chas*|d24cepyc02[5-9]|d24cepyc03[0-1])
#$ -N Dallas_FANOVA_Output

module purge
module load gsl
module load udunits
module load gdal
module load geos
module load R/4.2.1
module load gcc
Rscript Dallas_FANOVA_Postprocessing_PRCP_Jiachen.R

