FROM ubuntu:focal

RUN apt update && \
    apt upgrade -y && \
    apt install -y git pypy pypy-setuptools pypy-pytest

RUN git clone https://github.com/drtconway/zotmer.git && \
    cd zotmer && \
    pypy setup.py install

