
#create a list of the image files:
for i in $(seq 0 1 499); do echo "file './images/timestep_${i}.png'"; done > list.txt
#combine files as video
ffmpeg -r 20 -f concat -safe 0 -i list.txt -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -c:v libx264 -crf 15 -pix_fmt yuv420p compchem.mp4