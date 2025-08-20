#!/bin/csh
#$ -M mbonas@nd.edu
#$ -m ae
#$ -pe smp 24
#$ -q long
#$ -l h=(d12chas*|d24cepyc02[5-9]|d24cepyc03[0-1])
#$ -N DFW_PRCP_D01
#$ -t 1-1


module purge
module load gsl
module load udunits
module load gdal
module load geos
module load R
module load gcc
Rscript FANOVA_Dallas_CRC_PRCP.R


