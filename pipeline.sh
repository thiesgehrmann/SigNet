#!/bin/sh

./preprocessing.py
matlab < diffusion.m

./analysis.py

