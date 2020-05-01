
FROM alpine:latest as BUILD

RUN apk update && \
    apk add --no-cache build-base cmake

ENV SRC_PATH=/usr/src/ReDyMo-CPP

COPY . ${SRC_PATH}

WORKDIR ${SRC_PATH}/build

RUN cmake .. && \
    make -j2 && \
    make test

# Multi-stage build

FROM alpine:latest

ENV SRC_PATH=/usr/src/ReDyMo-CPP
ENV ORGANISM=TBrucei_TREU927
ENV APP_PATH=/opt/redymo

RUN addgroup -S redymo && adduser -S redymo -G redymo

RUN mkdir -p ${APP_PATH}; chown redymo: ${APP_PATH}

RUN apk update && \
    apk add --no-cache libstdc++ libgomp libgcc

USER redymo

WORKDIR ${APP_PATH}

COPY --from=BUILD ${SRC_PATH}/build/simulator ${APP_PATH}/

COPY --from=BUILD ${SRC_PATH}/data/database.sqlite ${APP_PATH}/data/

COPY --from=BUILD ${SRC_PATH}/data/MFA-Seq_${ORGANISM} ${APP_PATH}/data/MFA-Seq_${ORGANISM}

CMD ./simulator --cells 2 --organism 'Trypanosoma brucei brucei TREU927' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true --data-dir data
