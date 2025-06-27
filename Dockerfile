FROM rockylinux:9

RUN yum makecache --refresh
RUN yum -y install openmpi
RUN yum -y install make
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
RUN make all
RUN make mpi

#WORKDIR /Q/test/test_small
WORKDIR /Q/test/test_big

# Set the default command to perpetually running process
# ENTRYPOINT ["tail", "-f", "/dev/null"]

# Start running the tests right away
#ENTRYPOINT ["sh", "FEP_submit.sh"]
EXPOSE 5000