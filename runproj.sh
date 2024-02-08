
rm -rf output
rm -rf images
rm compchem.mp4

mkdir output
mkdir images

make re

./compchem

pvpython make_snapshots.py
sh make_movie.sh



