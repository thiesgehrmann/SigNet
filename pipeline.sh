#!/bin/sh

./preprocessing.py
matlab < diffusion.m


for beta in 0 0.0033 0.0067 0.0100 0.0133 0.0167 0.0200 0.0233 0.0267 0.0300; do
  sauto short --cpus-per-task 1 --mem 60000 "matlab -r \"beta=$beta;diffusion\""
done

./analysis.py

