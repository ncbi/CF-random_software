FROM nvidia/cuda:12.2.2-cudnn8-runtime-ubuntu22.04

LABEL maintainer="Pramesh Sharma"
LABEL description="CF-random: predicting alternative conformations and fold-switching proteins"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    bzip2 \
    ca-certificates \
    libgl1 \
    libglu1-mesa \
    libxi6 \
    libxrender1 \
    libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PYMOL_FLAGS="-c"
ENV CONDA_DIR=/opt/conda
ENV PATH="${CONDA_DIR}/bin:${PATH}"

RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p "${CONDA_DIR}" \
    && rm /tmp/miniforge.sh \
    && conda config --set solver libmamba \
    && conda clean -afy

WORKDIR /app
COPY environment.yml .

RUN conda env create -f environment.yml

ENV PATH="${CONDA_DIR}/envs/cf-random/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=cf-random

COPY . .
RUN pip install --no-cache-dir --retries 5 --timeout 100 -e .

RUN pip install --no-cache-dir --retries 5 --timeout 120 \
    "jax==0.4.28" "jaxlib==0.4.28+cuda12.cudnn89" \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

ENTRYPOINT ["cf-random"]
CMD ["--help"]