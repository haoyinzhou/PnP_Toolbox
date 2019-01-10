#!/bin/bash

ln -s ../data ./data 
matlab -nodisplay -nosoftwareopengl -r \
  "addpath(genpath('../code')); main_ordinary_no_outlier_npt; main_ordinary_no_outlier_sigma; main_ordinary_outliers_ransac; main_quasi_no_outlier_npt; main_quasi_no_outlier_sigma; main_quasi_outliers_ransac"