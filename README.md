# NLM
Non-local Means filter for 2D image de-noising
You can modify the main.cpp to add different types of noises like: Gaussian noise, impulse noise (NLM is not good at it) and Rician noise (has already added unbias).

Step:
mkdir build && cd build
cmake ..
make
./nlmClassic ../data/lena.jpg

