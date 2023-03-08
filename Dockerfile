FROM alpine:3.13.0 as BUILD

RUN apk update && \
    apk add --no-cache build-base cmake python3 py3-pip

ENV SRC_PATH=/usr/src/ReDyMo-CPP

COPY . ${SRC_PATH}

WORKDIR ${SRC_PATH}/build

RUN cmake .. && \
    make -j6 && \
    make test

# Multi-stage build
FROM alpine:3.13.0

ENV SRC_PATH=/usr/src/ReDyMo-CPP
ENV ORGANISM=TcruziCLBrenerEsmeraldo-like
ENV APP_PATH=/opt/redymo

RUN pip install numpy optuna

RUN addgroup -S redymo && adduser -S redymo -G redymo

RUN mkdir -p ${APP_PATH}; chown redymo: ${APP_PATH}

RUN apk update && \
    apk add --no-cache libstdc++ libgomp libgcc

USER redymo

WORKDIR ${APP_PATH}

COPY --from=BUILD ${SRC_PATH}/build/simulator ${APP_PATH}/

COPY --from=BUILD ${SRC_PATH}/data/database.sqlite ${APP_PATH}/data/

COPY --from=BUILD ${SRC_PATH}/data/MFA-Seq_${ORGANISM} ${APP_PATH}/data/MFA-Seq_${ORGANISM}

CMD ./simulator --cells 2 --organism '${ORGANISM}' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true --data-dir data
