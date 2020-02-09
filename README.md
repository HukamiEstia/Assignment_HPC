# HIGH PERFORMANCE COMPUTING
This repository contains the code produced for the High performance computers module at Cranfield University.

## INSTRUCTION
Build the docker image:

```docker build -t mpi-dev .```

start the container:

```docker run -it --rm -v $(pwd):/home mpi-dev:latest```

Inside the container:

```
cd home
make
mpirun -np <p> ./run <space> <time> <scheme>
```
