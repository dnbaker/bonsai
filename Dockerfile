FROM gcc:7

ADD . /code

WORKDIR /code

RUN make

RUN make unit
