FROM rootproject/root:6.30.02-alma9 AS mach3_build

LABEL maintainer="The MaCh3 Collaboration"
LABEL website="https://mach3-software.github.io/MaCh3/"
LABEL compiler="GNU 11.4.1"
LABEL root_version="v6.30.02"
LABEL org.opencontainers.image.description="Official MaCh3 container"

# install required packages
RUN git clone https://github.com/dpiparo/vdt.git /opt/vdt && \
    cd /opt/vdt && \
    cmake . && make install

RUN dnf install -y \
    diffutils

# install some common packages, can install more on demand
RUN dnf install -y \
    vim nano

# --- MaCh3 DUNE build ---

# Declare the build argument
ARG MACH3_DUNE_VERSION
ENV MACH3_DUNE_VERSION=${MACH3_DUNE_VERSION:-develop}

ENV MACH3_DUNE_WORK_DIR=/opt/MaCh3DUNE/
ENV MACH3_DUNE_INSTALL_DIR=${MACH3_DUNE_WORK_DIR}/build

ARG MACH3_DUNE_BUILD_ARGS

RUN git clone https://github.com/DUNE/MaCh3_DUNE.git -b ${MACH3_DUNE_VERSION} ${MACH3_DUNE_WORK_DIR} && \
    cd ${MACH3_DUNE_WORK_DIR} && \
    mkdir -p ${MACH3_DUNE_INSTALL_DIR} && \
    cd ${MACH3_DUNE_INSTALL_DIR} && \
    cmake ../ ${MACH3_DUNE_BUILD_ARGS} && \
    make -j && make install
RUN 

WORKDIR ${MACH3_DUNE_WORK_DIR}

ENTRYPOINT ["/bin/bash", "-c", "source build/bin/setup.MaCh3DUNE.sh && exec \"$@\"", "--"]
CMD ["/bin/bash"]