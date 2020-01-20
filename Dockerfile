FROM ubuntu:18.04

RUN apt update && apt install -y --no-install-recommends \
    build-essential \
    hwloc-nox \
    libc6 \
    libcr0 \
    libhwloc5 \
    libmpich12 \
    libmpich-dev \
    mpich \
    openssh-client \
    openssh-server \