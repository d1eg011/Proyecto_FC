#!/bin/bash

./a.out
# Create a temporary directory for intermediate files
mkdir -p temp_images

# Ensure the ordering of the PDFs by extracting the numeric part and sorting numerically
count=1
for pdf in $(ls lattice_*.pdf | sort -t_ -k2,2n); do
    # Convert the PDF to a PNG image
    pdftoppm -png -singlefile "$pdf" "temp_images/img$count"
    count=$((count + 1))
done

# Create a video from the PNG files
ffmpeg -framerate 10 -i temp_images/img%d.png -c:v libx264 -r 30 -pix_fmt yuv420p Simulation_Ising.mp4

# Clean up temporary images
rm -r temp_images *lattice* *.out

