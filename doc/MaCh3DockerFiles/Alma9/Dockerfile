#To run use: docker build --secret id=tokens,src=tokens.txt -t mach3dune .
FROM picker24/root_v6_26_10:alma9 AS mach3dune_build

# Declare the build argument
ARG MACH3_DUNE_VERSION

ENV MACH3_DUNE_VERSION=${MACH3_DUNE_VERSION:-develop}

ENV MACH3_DUNE_WORK_DIR=/opt/MaCh3DUNE/
ENV MACH3_DUNE_BUILD_DIR=${MACH3_DUNE_WORK_DIR}/build
ENV MACH3_DUNE_INSTALL_DIR=${MACH3_DUNE_WORK_DIR}/install


ARG MACH3_DUNE_BUILD_ARGS

RUN git clone https://github.com/DUNE/MaCh3_DUNE.git ${MACH3_DUNE_WORK_DIR}

WORKDIR ${MACH3_DUNE_WORK_DIR}

RUN git checkout ${MACH3_DUNE_VERSION}

# Make a python environment
RUN pip install virtualenv
RUN virtualenv .env
RUN source .env/bin/activate
RUN pip install -r requirements.txt

RUN mkdir -p ${MACH3_DUNE_BUILD_DIR}
WORKDIR ${MACH3_DUNE_BUILD_DIR}

RUN cmake -DCMAKE_INSTALL_PREFIX=${MACH3_DUNE_INSTALL_DIR} ${MACH3_DUNE_BUILD_ARGS} ..
RUN cmake --build . --target install -j8

WORKDIR ${MACH3_DUNE_WORK_DIR}
