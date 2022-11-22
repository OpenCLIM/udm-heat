FROM python:3.11

RUN mkdir /src
RUN mkdir /data

WORKDIR /src

RUN apt-get -y update

RUN apt-get -y install libgdal-dev gdal-bin

COPY main.py /src
#COPY data /data

ENTRYPOINT python main.py
