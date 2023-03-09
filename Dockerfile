FROM ubuntu:20.04 as BUILD

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake

ENV SRC_PATH=/usr/src/ReDyMo-CPP

COPY . ${SRC_PATH}

WORKDIR ${SRC_PATH}/build

RUN cmake .. -DCMAKE_BUILD_TYPE=Performance -DBUILD_TESTING=ON -DCOVERAGE=OFF -DBUILD_GPGPU=OFF && \
    make -j15 && \
    make test

# Multi-stage build
FROM ubuntu:20.04

ENV SRC_PATH=/usr/src/ReDyMo-CPP
ENV ORGANISM=TcruziCLBrenerEsmeraldo-like
ENV APP_PATH=/opt/redymo

# RUN addgroup -S redymo && adduser -S redymo -G redymo

# RUN mkdir -p ${APP_PATH}; chown redymo: ${APP_PATH}

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake python3 python3-pip

RUN pip install numpy optuna pandas

# USER redymo

WORKDIR ${APP_PATH}

COPY --from=BUILD ${SRC_PATH}/build/simulator ${APP_PATH}/

COPY --from=BUILD ${SRC_PATH}/data/database.sqlite ${APP_PATH}/data/

COPY --from=BUILD ${SRC_PATH}/data/MFA-Seq_${ORGANISM} ${APP_PATH}/data/MFA-Seq_${ORGANISM}

COPY --from=BUILD ${SRC_PATH}/script ${APP_PATH}/script

CMD ./simulator --cells 2 --organism '${ORGANISM}' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true --data-dir data
