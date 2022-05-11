FROM fedora

# Install some basic utilities
# RUN apt-get update && apt-get install -y \
#     curl \
#     ca-certificates \
#     sudo \
#     git \
#     bzip2 \
#     libx11-6 \
#     gcc-7 g++-7 \
#     gfortran-7 \
#  && rm -rf /var/lib/apt/lists/*

RUN yum install -y gcc-gfortran

# Create a non-root user and switch to it
# RUN adduser --disabled-password --gecos '' --shell /bin/bash user
# RUN echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/90-user

# Create a working directory

RUN mkdir /Q
WORKDIR /Q

# RUN chown -R user /home/user/Q

USER root

# All users can use /home/user as their home directory
ENV HOME=/Q
RUN chmod 777 /Q

COPY /bin/q6 ./bin
COPY /test/q6 ./test

# Set the default command to perpetually running process
ENTRYPOINT ["tail", "-f", "/dev/null"]