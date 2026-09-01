FROM rootproject/root:6.32.02-ubuntu22.04 AS mach3_dune_build

LABEL maintainer="The MaCh3 Collaboration"
LABEL website="https://dune.github.io/MaCh3_DUNE/"
LABEL compiler="GNU 11.4.0"
LABEL root_version="v6.32.02"
LABEL org.opencontainers.image.description="Official MaCh3 DUNE container"

# install required packages
RUN apt update && apt upgrade -y
RUN apt install -y --no-install-recommends nlohmann-json3-dev

# install some common packages, can install more on demand
RUN apt install -y --no-install-recommends \
   vim nano    

# --- MaCh3 DUNE build ---

# Declare the build argument
ARG MACH3_DUNE_VERSION
ENV MACH3_DUNE_VERSION=${MACH3_DUNE_VERSION:-develop}

ENV MACH3_DUNE_WORK_DIR=/opt/MaCh3DUNE/
ENV MACH3_DUNE_INSTALL_DIR=${MACH3_DUNE_WORK_DIR}/build

ARG MACH3_DUNE_BUILD_ARGS

RUN git clone https://github.com/DUNE/MaCh3_DUNE.git ${MACH3_DUNE_WORK_DIR} && \
    cd ${MACH3_DUNE_WORK_DIR} && \
    git fetch origin ${MACH3_DUNE_VERSION} && \
    git checkout ${MACH3_DUNE_VERSION} && \
    mkdir -p ${MACH3_DUNE_INSTALL_DIR} && \
    cd ${MACH3_DUNE_INSTALL_DIR} && \
    cmake ../ ${MACH3_DUNE_BUILD_ARGS} && \
    make -j && make install

WORKDIR ${MACH3_DUNE_WORK_DIR}

ENTRYPOINT ["/bin/bash", "-c", "source build/bin/setup.MaCh3DUNE.sh && exec \"$@\"", "--"]
CMD ["/bin/bash"]