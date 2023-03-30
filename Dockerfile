# Multi-stage build
FROM ubuntu:20.04 as EXECUTER

ENV SRC_PATH=/usr/src/ReDyMo-CPP
ENV ORGANISM=TcruziCLBrenerEsmeraldo-like
ENV APP_PATH=/opt/redymo

# RUN addgroup -S redymo && adduser -S redymo -G redymo

# RUN mkdir -p ${APP_PATH}; chown redymo: ${APP_PATH}

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake python3 python3-dev python3-pip vim sqlite3 default-libmysqlclient-dev

RUN pip install numpy optuna pandas mysqlclient

FROM ubuntu:20.04 as BUILDER

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential cmake

FROM BUILDER as COMPILER

ENV SRC_PATH=/usr/src/ReDyMo-CPP

COPY . ${SRC_PATH}

WORKDIR ${SRC_PATH}/build

RUN cmake .. -DCMAKE_BUILD_TYPE=Performance -DBUILD_TESTING=ON -DCOVERAGE=OFF -DBUILD_GPGPU=OFF && \
    make -j15 && \
    make test

FROM EXECUTER

# USER redymo

WORKDIR ${APP_PATH}

COPY --from=COMPILER ${SRC_PATH}/build/simulator ${APP_PATH}/

COPY --from=COMPILER ${SRC_PATH}/data/database.sqlite ${APP_PATH}/data/

COPY --from=COMPILER ${SRC_PATH}/data/MFA-Seq_${ORGANISM} ${APP_PATH}/data/MFA-Seq_${ORGANISM}

COPY --from=COMPILER ${SRC_PATH}/script ${APP_PATH}/script

VOLUME ${APP_PATH}/train-db

WORKDIR /opt/redymo/script

CMD nice -n 20 nohup python3 train-model-hyperparameters.py
