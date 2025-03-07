#To run use: docker build --secret id=tokens,src=tokens.txt -t mach3dune .
FROM rootproject/root:6.32.02-ubuntu22.04 AS mach3_dune_build

RUN apt update && apt upgrade -y
RUN apt-get install -y nlohmann-json3-dev
#MISC_SW
RUN apt install -y --no-install-recommends \
    vim less nano gdb csh tcsh ed quota python3 python3-dev python3-pip \
    cvs procmail ca-certificates cmake ninja-build


# Declare the build argument
ARG MACH3_DUNE_VERSION
ENV MACH3_DUNE_VERSION=${MACH3_DUNE_VERSION:-develop}

ENV MACH3_DUNE_WORK_DIR=/opt/MaCh3DUNE/
ENV MACH3_DUNE_INSTALL_DIR=${MACH3_DUNE_WORK_DIR}/build


RUN --mount=type=ssh git clone https://github.com/DUNE/MaCh3_DUNE.git ${MACH3_DUNE_WORK_DIR}

WORKDIR ${MACH3_DUNE_WORK_DIR}
RUN git checkout ${MACH3_DUNE_VERSION}

RUN mkdir -p ${MACH3_DUNE_INSTALL_DIR}
WORKDIR ${MACH3_DUNE_INSTALL_DIR}

RUN cmake ../

RUN make -j && make install