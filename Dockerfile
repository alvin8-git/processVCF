# =============================================================================
# Dockerfile for processVCF Pipeline
# =============================================================================
# Lightweight Ubuntu-based container with all required annotation tools
#
# Build: docker build -t processvcf .
# Run:   docker run -v /path/to/Databases:/home/user/Databases \
#                   -v /path/to/data:/data \
#                   processvcf
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="Alvin Ng"
LABEL description="VCF Processing Pipeline for TMSP and CEBPA panels"
LABEL version="2.1"

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Create non-root user
ARG USERNAME=user
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME

# =============================================================================
# SYSTEM PACKAGES
# =============================================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
    # Core utilities
    bash \
    coreutils \
    gawk \
    sed \
    grep \
    curl \
    wget \
    git \
    ca-certificates \
    # Bioinformatics tools
    bcftools \
    tabix \
    vcftools \
    # Programming languages
    perl \
    python3 \
    python3-pip \
    openjdk-11-jre-headless \
    # Parallel processing
    parallel \
    # Compression
    gzip \
    bzip2 \
    xz-utils \
    unzip \
    # Required for Perl modules and compilation
    build-essential \
    cpanminus \
    zlib1g-dev \
    python3-dev \
    # Required for VEP
    libdbi-perl \
    libdbd-mysql-perl \
    libwww-perl \
    libjson-perl \
    libarchive-zip-perl \
    # Rename utility
    rename \
    && rm -rf /var/lib/apt/lists/*

# =============================================================================
# PERL MODULES
# =============================================================================

RUN cpanm --notest Excel::Writer::XLSX

# =============================================================================
# PYTHON PACKAGES (TransVar)
# =============================================================================

RUN pip3 install --no-cache-dir transvar

# Configure TransVar with hg19 reference
RUN mkdir -p /home/$USERNAME/.transvar.cfg.d

# =============================================================================
# SOFTWARE DIRECTORY SETUP
# =============================================================================

# Create directories for software and databases
RUN mkdir -p /home/$USERNAME/Software \
             /home/$USERNAME/Databases \
             /home/$USERNAME/Scripts

# =============================================================================
# ANNOVAR (copy from local or download)
# Uncomment the COPY line if you have ANNOVAR locally
# =============================================================================

# ANNOVAR requires registration - copy from local installation
COPY --chown=$USERNAME:$USERNAME annovar/ /home/$USERNAME/Software/annovar/

# =============================================================================
# snpEff (copy from local)
# Note: snpEff database should be mounted at runtime
# Mount path: /home/user/Databases/snpEff
# =============================================================================

COPY --chown=$USERNAME:$USERNAME snpEff/ /home/$USERNAME/Software/snpEff/

# Create symlink for snpEff data directory (will be mounted at runtime)
RUN mkdir -p /home/$USERNAME/Databases/snpEff && \
    ln -sf /home/$USERNAME/Databases/snpEff /home/$USERNAME/Software/snpEff/data

# =============================================================================
# CancerVar (copy from local)
# =============================================================================

COPY --chown=$USERNAME:$USERNAME CancerVar/ /home/$USERNAME/Software/CancerVar/

# =============================================================================
# VEP (Ensembl Variant Effect Predictor)
# =============================================================================
# VEP installation is optional - can be mounted from host or installed separately
# Mount path for VEP: /home/user/Software/ensembl-vep
# Mount path for VEP cache: /home/user/Databases/vep

# Create VEP directories for mounting
RUN mkdir -p /home/$USERNAME/Software/ensembl-vep \
             /home/$USERNAME/Databases/vep

# =============================================================================
# PIPELINE SCRIPTS
# =============================================================================

COPY --chown=$USERNAME:$USERNAME processVCF.sh /home/$USERNAME/Scripts/
COPY --chown=$USERNAME:$USERNAME mergeVCFannotation-optimized.sh /home/$USERNAME/Scripts/

RUN chmod +x /home/$USERNAME/Scripts/*.sh

# =============================================================================
# ENVIRONMENT VARIABLES
# =============================================================================

ENV HOME=/home/$USERNAME
ENV PATH="/home/$USERNAME/Software/annovar:/home/$USERNAME/Software/snpEff:/home/$USERNAME/Software/ensembl-vep:/home/$USERNAME/Scripts:$PATH"
ENV PERL5LIB="/home/$USERNAME/Software/ensembl-vep/modules"

# =============================================================================
# WORKING DIRECTORY AND USER
# =============================================================================

WORKDIR /data
USER $USERNAME

# Default command
CMD ["bash"]
