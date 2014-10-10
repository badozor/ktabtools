#!/bin/sh 
R CMD check ktabtools --library="/export/scratch/R/library"
R CMD build ktabtools --library="/export/scratch/R/library"
R CMD build -binary ktabtools --library="/export/scratch/R/library"
R CMD Rd2pdf --pdf ktabtools --library="/export/scratch/R/library"
