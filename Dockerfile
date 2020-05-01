
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

RUN apk update && \
    apk add --no-cache libstdc++ libgomp libgcc

RUN useradd -ms /bin/bash ReDyMo

USER ReDyMo

ENV SRC_PATH=/usr/src/ReDyMo-CPP

ENV ORGANISM=TBrucei_TREU927

COPY --from=BUILD ${SRC_PATH}/build/simulator /opt/redymo/

COPY --from=BUILD ${SRC_PATH}/data/simulation.sqlite /opt/redymo/data/

COPY --from=BUILD ${SRC_PATH}/data/MFA-Seq_${ORGANISM} /opt/redymo/data/MFA-Seq_${ORGANISM}

CMD ls /opt
