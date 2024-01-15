FROM mambaorg/micromamba:1.5.6
# using base from https://github.com/mamba-org/micromamba-docker

ARG BAKTA_VERSION=1.9.1
ARG BWA_VERSION=0.7.17
ARG LIBDEFLATE_VERSION=1.18
ARG SAMTOOLS_VERSION=1.19

ENV PATH="/opt/conda/bin:${PATH}"

USER root

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install -y --no-install-recommends wget tar

RUN apt-get update -y && apt-get upgrade -y && \
    apt-get install -y python-is-python3

USER $MAMBA_USER

RUN wget -O /tmp/bakta-${BAKTA_VERSION}.tar.gz https://github.com/oschwengers/bakta/archive/refs/tags/v${BAKTA_VERSION}.tar.gz \
    && tar -xzf /tmp/bakta-${BAKTA_VERSION}.tar.gz \
    && rm /tmp/bakta-${BAKTA_VERSION}.tar.gz \
    && chown $MAMBA_USER:$MAMBA_USER /tmp/bakta-${BAKTA_VERSION}/

RUN echo -e '\n  - deepsig>=1.2.5' >> /tmp/bakta-${BAKTA_VERSION}/environment.yml

RUN micromamba install -y -n base -f /tmp/bakta-${BAKTA_VERSION}/environment.yml && \
    micromamba install -c bioconda -y -n base bwa=${BWA_VERSION} && \
    micromamba install -c conda-forge -y -n base libdeflate=${LIBDEFLATE_VERSION} && \
    micromamba install -c bioconda -y -n base samtools=${SAMTOOLS_VERSION} && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN python3 -m pip install -v --no-cache /tmp/bakta-${BAKTA_VERSION}/

COPY requirements.txt .

RUN pip install -r requirements.txt