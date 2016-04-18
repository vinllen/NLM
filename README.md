# NLM
## Non-local Means filter for 2D image de-noising
You can modify the main.cpp to add different level of noise and types of noises like: Gaussian noise, impulse noise (NLM is not good at it) and Rician noise (has already added unbias). <br \>

## Step:
mkdir build && cd build <br \>
cmake .. <br \>
make <br \>
./nlmClassic ../data/lena.jpg <br \>

