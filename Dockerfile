FROM ubuntu:latest
RUN apt-get update
RUN apt-get install python3 trimmomatic fastqc
WORKDIR Software
RUN wget https://sourceforge.net/projects/bbmap/files/BBMap_39.01.tar.gz/download
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip
