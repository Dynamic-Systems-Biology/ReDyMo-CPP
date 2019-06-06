FROM rikorose/gcc-cmake:latest

ENV SRC_PATH=/usr/src/ReDyMo-CPP

WORKDIR ${SRC_PATH}

COPY . .

WORKDIR ./build

RUN cmake .. && \
    make && \
    make test

CMD ./simulator
