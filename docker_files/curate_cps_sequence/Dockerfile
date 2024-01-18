FROM ubuntu:23.04

RUN apt-get update -qq -y && apt-get upgrade -qq -y && \
    apt-get install -y  \
    bedtools \
    python-is-python3 \
    python3-pip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


COPY requirements.txt .

RUN pip3 install -r requirements.txt --break-system-packages && \
    rm requirements.txt
