FROM spiralaico/python38_pytorch_base

RUN pip install torch
RUN pip install networkx

RUN pip install deepsvp 

RUN git clone https://github.com/lgmgeo/AnnotSV.git --branch v2.3
RUN cd AnnotSV && make PREFIX=. install
RUN cd AnnotSV && make DESTDIR= PREFIX=. install-human-annotation

ENV ANNOTSV=/AnnotSV
