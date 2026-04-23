FROM docker.io/bioconductor/bioconductor_docker:RELEASE_3_22

# System deps:
#   libhdf5-dev  : h5ad / h5seurat / TENxPBMCData
#   python3-*    : scanpy side
#   rust         : scx (only distributed via cargo)
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        python3 \
        python3-pip \
        python3-venv \
        pkg-config \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Minimal rust toolchain for `cargo install scx-cli`.
ENV CARGO_HOME=/usr/local/cargo RUSTUP_HOME=/usr/local/rustup
ENV PATH=${CARGO_HOME}/bin:${PATH}
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --profile minimal --default-toolchain stable \
 && cargo install --locked scx-cli

RUN cargo install --features gpu denet
# RUN python3 -m pip install --no-cache-dir --break-system-packages denet

# R packages: Bioconductor + Seurat + anndataR + BPCells.
RUN Rscript -e " \
    BiocManager::install(c( \
        'SingleCellExperiment', \
        'DelayedArray', \
        'TENxPBMCData', \
        'TENxBrainData', \
        'anndataR' \
    ), ask = FALSE, update = FALSE); \
    install.packages(c('argparse', 'Seurat', 'remotes'), repos = 'https://cloud.r-project.org'); \
    remotes::install_github('bnprks/BPCells/r', upgrade = 'never') \
    "

# Python: scanpy + anndata.
# --ignore-installed avoids pip trying to uninstall debian-managed packages
# (e.g. packaging, threadpoolctl) which lack RECORD files.
RUN python3 -m pip install --no-cache-dir --break-system-packages --ignore-installed \
        scanpy \
        anndata \
        h5py

# Prevent BLAS/libgomp from claiming all cores at thread-pool init time.
# Per-call pragmas still override these where real parallelism is intended.
ENV OMP_NUM_THREADS=1 \
    OMP_THREAD_LIMIT=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# omnibench-logger: R and Python both live in btraven00/omnibench-logger.
RUN Rscript -e "remotes::install_github('btraven00/omnibench-logger', subdir = 'omnibench.logger', upgrade = 'never')"
RUN python3 -m pip install --no-cache-dir --break-system-packages --ignore-installed \
        pip setuptools wheel \
 && python3 -m pip install --no-cache-dir --break-system-packages \
        git+https://github.com/btraven00/omnibench-logger.git
