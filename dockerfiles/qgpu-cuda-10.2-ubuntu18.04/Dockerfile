FROM nvidia/cuda:10.2-base-ubuntu18.04

# Install some basic utilities
RUN apt-get update && apt-get install -y \
    curl \
    ca-certificates \
    sudo \
    git \
    bzip2 \
    libx11-6 \
    gcc-7 g++-7 \
    gfortran-7 \
 && rm -rf /var/lib/apt/lists/*

# Create a non-root user and switch to it
RUN adduser --disabled-password --gecos '' --shell /bin/bash user
RUN echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-user

# Create a working directory
RUN mkdir /home/user/Q
WORKDIR /home/user/Q

RUN chown -R user /home/user/Q

USER user

# All users can use /home/user as their home directory
ENV HOME=/home/user
RUN chmod 777 /home/user

# Install Miniconda and Python 3.8
ENV CONDA_AUTO_UPDATE_CONDA=false
ENV PATH=/home/user/miniconda/bin:$PATH
RUN curl -sLo ~/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-py38_4.8.3-Linux-x86_64.sh \
 && chmod +x ~/miniconda.sh \
 && ~/miniconda.sh -b -p ~/miniconda \
 && rm ~/miniconda.sh \
 && conda install -y python==3.8.3 \
 && conda clean -ya

COPY /bin ./bin
COPY /src ./src
COPY /share ./share
COPY /test ./test
COPY /env ./env
COPY ./setup.py ./setup.py
COPY ./pyproject.toml ./pyproject.toml

RUN pip install .
RUN pip install numpy

# Set the default command to perpetually running process
ENTRYPOINT ["tail", "-f", "/dev/null"]