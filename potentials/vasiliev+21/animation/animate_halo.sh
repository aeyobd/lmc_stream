
set -xe 

animate_dm.jl projected_halo.hdf5 --scalebar 100

cd figures/dm_animation
ffmpeg -framerate 30 -i frame_%d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4
