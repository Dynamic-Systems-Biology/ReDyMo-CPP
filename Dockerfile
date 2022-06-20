
FROM alpine:3.13.0 as BUILD

RUN apk update && \
    apk add --no-cache build-base cmake bash pcc-libs-dev

ENV SRC_PATH=/usr/src/ReDyMo-CPP

COPY . ${SRC_PATH}

WORKDIR ${SRC_PATH}/build

RUN cmake .. -DGPU_ENABLED=false -DCMAKE_BUILD_TYPE=Performance -DCMAKE_VERBOSE_MAKEFILE=OFF && \
    make -j8 -s && \
    make test

# Multi-stage build
FROM alpine:3.13.0

ENV SRC_PATH=/usr/src/ReDyMo-CPP
ENV ORGANISM='Trypanosoma brucei brucei TREU927'
ENV APP_PATH=/opt/redymo

RUN addgroup -S redymo && adduser -S redymo -G redymo

RUN mkdir -p ${APP_PATH}; chown redymo: ${APP_PATH}

RUN apk update && \
    apk add --no-cache libstdc++ libgomp libgcc ruby

RUN gem install ruby-progressbar optparse

USER redymo

WORKDIR ${APP_PATH}

COPY --from=BUILD --chown=redymo:redymo ${SRC_PATH} ${APP_PATH}/

#COPY --from=BUILD ${SRC_PATH}/data/database.sqlite ${APP_PATH}/data/

#COPY --from=BUILD ${SRC_PATH}/script ${APP_PATH}/script/

#COPY --from=BUILD ${SRC_PATH}/data/MFA-Seq_${ORGANISM} ${APP_PATH}/data/MFA-Seq_${ORGANISM}

CMD ./simulator --cells 2 --organism 'Trypanosoma brucei brucei TREU927' --resources 10 --speed 65 --period 150 --timeout 1000000 --dormant true --data-dir data
