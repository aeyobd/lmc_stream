
set -xe 

animate_dm.jl projected_lmc.hdf5 --scalebar 100 -o lmc_animation

cd lmc_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
