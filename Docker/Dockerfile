FROM nfcore/base:1.12.1
LABEL authors="Carlos Ruiz-Arenas" \
      description="Docker image containing all software requirements for the Obesity Analysis"

# Instal basic Ubuntu programs
RUN apt-get update --fix-missing && \
  apt-get install -q -y build-essential

# Install mamba to manage conda
RUN conda install mamba -c conda-forge
# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/obesity_analysis/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name obesity_analysis > obesity_analysis.yml

# Install PAINTOR
RUN cd /opt/conda/ && git clone https://github.com/gkichaev/PAINTOR_V3.0.git && \
  cd PAINTOR_V3.0 && bash install.sh

# Install Bioconductor
RUN R -e 'BiocManager::install(c("Biostrings", "snpStats", "VariantAnnotation", "GenomicRanges", "BiocParallel", "SummarizedExperiment",\
                                  "SNPassoc"))'

# Install scoreInvHap
RUN R -e 'devtools::install_github("isglobal-brge/scoreInvHap")'

# Install libcurl for liftOver
RUN apt-get install -q -y libcurl4
