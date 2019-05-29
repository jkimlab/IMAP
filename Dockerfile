############################################################
# Dockerfile for Integrative Meta-Assembly Pipeline(IMAP)
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu 16.04
#FROM ubuntu:16.04
FROM ubuntu@sha256:e4a134999bea4abb4a27bc437e6118fdddfb172e1b9d683129b74d254af51675

# File Author / Maintainer
MAINTAINER Ho Yong Lee

# Update the repository sources list
RUN apt-get update && apt-get install -y \
    git \
    g++ \
    build-essential \
    pkg-config \
    libgd-dev \
    libncurses-dev \
    libghc-bzlib-dev \
    libboost-all-dev \
    python \ 
    wget \
    perl \
    cpanminus \
    openjdk-8-jdk \
    openjdk-8-jre

RUN cpanm \
    YAML \
    Switch \
    Parallel::ForkManager \
    ExtUtils::PkgConfig \
    GD \
    XML::Parser \
    XML::Parser::PerlSAX \
    XML::DOM \
    XML::DOM::XPath \
    XML::Twig \
    Bio::TreeIO

RUN cd home && git clone https://github.com/jkimlab/IMAP.git && \
    cd IMAP && perl build.pl --install 
