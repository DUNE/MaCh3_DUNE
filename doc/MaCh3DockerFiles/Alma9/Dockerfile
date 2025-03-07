#To run use: docker build --secret id=tokens,src=tokens.txt -t mach3dune .
FROM picker24/root_v6_26_10:alma9 AS mach3dune_build

# Declare the build argument
ARG MACH3_DUNE_VERSION
ARG MACH3_CORE_VERSION

ENV MACH3_DUNE_VERSION=${MACH3_DUNE_VERSION:-develop}

ENV MACH3_DUNE_WORK_DIR=/opt/MaCh3DUNE/
ENV MACH3_DUNE_INSTALL_DIR=${MACH3_DUNE_WORK_DIR}/build

ARG MACH3_DUNE_BUILD_ARGS

RUN git clone https://github.com/DUNE/MaCh3_DUNE.git ${MACH3_DUNE_WORK_DIR}

WORKDIR ${MACH3_DUNE_WORK_DIR}
RUN git checkout ${MACH3_DUNE_VERSION}

RUN mkdir -p ${MACH3_DUNE_INSTALL_DIR}
WORKDIR ${MACH3_DUNE_INSTALL_DIR}

RUN if [ -n "$MACH3_CORE_VERSION" ]; then \
      cmake ../ -DMaCh3_Branch=${MACH3_CORE_VERSION}; \
    else \
      cmake ../; \
    fi

RUN make -j && make install