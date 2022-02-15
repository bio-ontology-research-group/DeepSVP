FROM python:3.8

RUN apt update
RUN apt install -y bedtools bcftools

RUN pip install torch networkx

RUN git clone https://github.com/lgmgeo/AnnotSV.git --branch v2.3
WORKDIR /AnnotSV
RUN make PREFIX=. install

ENV ANNOTSV=/AnnotSV

RUN pip install deepsvp