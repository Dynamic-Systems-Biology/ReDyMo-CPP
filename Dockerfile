FROM rikorose/gcc-cmake:latest
COPY . /usr/src/ReDyMo-CPP
RUN mkdir build
WORKDIR /usr/src/ReDyMo-CPP/build
RUN cmake ..
RUN make
RUN make test
CMD ["./simulator"]
