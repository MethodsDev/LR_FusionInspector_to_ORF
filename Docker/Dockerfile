FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends \
    build-essential \
    ca-certificates \
    cpanminus \
    curl \
    diamond-aligner \
    git \
    hmmer \
    libbz2-dev \
    libdb-dev \
    liblzma-dev \
    libncurses-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncbi-blast+ \
    perl \
    python3 \
    python3-pip \
    python3-setuptools \
    samtools \
    tabix \
    unzip \
    wget \
    zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/bin/python3 /usr/bin/python && \
    cpanm --notest DB_File URI::Escape

RUN python3 -m pip install --break-system-packages \
    intervaltree \
    lmdb \
    networkx \
    psutil \
    pysam \
    tqdm

ENV APP_HOME=/opt/LR_FusionInspector_to_ORF
ENV LRAA_DIR=${APP_HOME}/LRAA
ENV TRANSDECODER_DIR=${APP_HOME}/TransDecoder

WORKDIR ${APP_HOME}

COPY LR_FusionInspector_to_ORF.py ${APP_HOME}/
COPY LRAA ${LRAA_DIR}
COPY TransDecoder ${TRANSDECODER_DIR}

RUN chmod +x ${APP_HOME}/LR_FusionInspector_to_ORF.py \
    ${LRAA_DIR}/LRAA \
    ${TRANSDECODER_DIR}/TransDecoder && \
    make -C ${TRANSDECODER_DIR}

ENV PATH="${APP_HOME}:${LRAA_DIR}:${LRAA_DIR}/util:${LRAA_DIR}/util/misc:${TRANSDECODER_DIR}:${PATH}"

ENTRYPOINT ["python3", "/opt/LR_FusionInspector_to_ORF/LR_FusionInspector_to_ORF.py"]
