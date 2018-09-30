# Base setup is an anaconda python2.7 environment
FROM continuumio/anaconda

# Run code to do the following:
# 1. Set up a new user directory with python and git subdirectories
# 2. Set up a directory for the cdfassnacht matplotlibrc
# 3. Add a user called tmpuser
# 
RUN mkdir /opt/notebooks && \
    useradd -ms /bin/bash tmpuser && \
    mkdir -p /home/tmpuser/git && \
    mkdir -p /home/tmpuser/python && \
    /opt/conda/bin/conda update jupyter

# Put the matplotlibrc file that is included with this image into the
#  /home/tmpuser/.matplotlib directory
COPY matplotlibrc /home/tmpuser/

# Set up the user of the container to be the new tmpuser
USER tmpuser
WORKDIR /home/tmpuser

# Run a command
#CMD ipython
CMD ["/opt/conda/bin/jupyter", "notebook", "--notebook-dir=/home/tmpuser", "--ip='*'", "--port=8888", "--no-browser"]
