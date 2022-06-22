FROM rockylinux

RUN yum makecache --refresh
# RUN yum -y install openmpi
RUN yum -y install make
COPY /dockerfiles/q-ifort/oneAPI.repo /tmp/oneAPI.repo
RUN mv /tmp/oneAPI.repo /etc/yum.repos.d
RUN yum install -y intel-hpckit
ENV PATH=/usr/lib64/openmpi/bin/:$PATH

# Create a working directory

RUN mkdir /Q
WORKDIR /Q

USER root

# All users can use /home/user as their home directory
ENV HOME=/Q
RUN chmod 777 /Q

# compile Q from source
COPY /src/q6 ./src/q6
COPY /test/q6 ./test

RUN mkdir bin

WORKDIR /Q/src/q6
RUN make all COMP=ifort
RUN make mpi COMP=ifort

WORKDIR /Q

# Set the default command to perpetually running process
ENTRYPOINT ["tail", "-f", "/dev/null"]