#!/bin/bash

echo "-------------------------"
echo "ACAQ visualizer for Linux"
echo "-------------------------"

command -v R >/dev/null 2>&1 || { echo >&2 "Please install R (Debian/Ubuntu: sudo apt install R-base)"; exit 1; }

cd ACAQ_Shiny

# Packrat initialization part 1
echo "[1/3] Bootstrapping packrat package management system ..."
R --vanilla -f packrat/init.R --args --bootstrap-packrat

# Force packrat to REALLY restore all libraries (because sometimes it does NOT)
echo "[2/3] Forcing packrat to restore all libraries ..."
echo "
source('packrat/init.R')
packrat::restore()
" | R --vanilla

# Start ACAQ Visualizer
echo "
source('packrat/init.R')
shiny::runApp(launch.browser=T)
" | R --vanilla
