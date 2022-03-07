FROM continuumio/miniconda3:latest

LABEL maintainer="Sebastian Uhrig <s.uhrig@dkfz.de>"

SHELL ["/bin/bash", "-c"]

# install software using conda
COPY task-environment.yml ./
RUN conda init bash && \
    conda env create -n main -f task-environment.yml && \
    source activate main && \
    conda clean --all -f -y

# ps is needed for collecting runtime information from the container
RUN apt update && \
    apt-get install -y procps && \
    rm -rf /var/lib/apt/lists/*

# For login Bash /etc/profile and ~/.profile is sourced. /etc/profile sources /etc/bash.bashrc.
# For non-login, interactive Bash /etc/bash.bashrc is sourced directly.
# For non-login, non-interactive Bash. We set BASH_ENV/ENV to /etc/bash.bashrc
# NOTE: For unknown reasons /.bashrc could not be used, because when using
#       `-u $(id -u):$(id -g)` as docker run parameter, the file was absent
#       (but not when just starting the container as root). Therefore
#       /etc/bash.bashrc is used.
# NOTE: Conda should be fully available in non-login, interactive shell. Conda itself creates
#       /etc/profile.d/conda.sh. The code that `conda init bash` writes to ~/.bashrc is moved
#       to /etc/bash.bashrc and reads the /etc/profile.d/conda.sh.
ENV BASH_ENV /etc/container.bashrc
ENV ENV /etc/container.bashrc

RUN grep "managed by 'conda init'" -A 100 ~/.bashrc >> /etc/container.bashrc && \
    rm ~/.bashrc && \
    echo -e '\
set +u\n\
source activate main\n\
set -u\n' >> /etc/container.bashrc && \
    echo "source /etc/profile" > ~/.profile && \
    cp ~/.profile /.profile && \
    echo "source /etc/container.bashrc" >> /etc/bash.bashrc

ENTRYPOINT ["bash", "-i", "-c"]
