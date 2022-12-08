FROM python:3.11

RUN mkdir /src
RUN mkdir /data

WORKDIR /src

RUN apt-get -y update

RUN apt-get -y install libgdal-dev gdal-bin

RUN pip3 install geopandas rasterio

COPY main.py /src

ENTRYPOINT python main.py
