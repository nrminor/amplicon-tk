FROM ubuntu:20.04

# Set working directory
WORKDIR /scratch

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York
ENV HOME /home

# Set default command to be the bash shell
ENTRYPOINT ["bash"]

# Install a few ubuntu dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    curl \
    wget \
    make \
    gcc \
    cmake \
    git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir /dependencies && \
    dpkg -l > /dependencies/apt-get.lock

# install just
RUN mkdir -p ~/bin && \
    curl --proto '=https' --tlsv1.2 -sSf https://just.systems/install.sh | bash -s -- --to ~/bin
ENV PATH "$PATH:$HOME/bin"

# run just environment recipe
RUN just env
ENV PATH "$PATH:/.venv/bin"
