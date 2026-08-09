FROM condaforge/miniforge3:24.11.0-0

WORKDIR /opt/coalescent-based-dating

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /opt/coalescent-based-dating/environment.yml
RUN conda env create -f /opt/coalescent-based-dating/environment.yml && conda clean -afy

COPY . /opt/coalescent-based-dating
RUN /bin/bash -lc "source /opt/conda/etc/profile.d/conda.sh && conda activate cbdating && bash scripts/bootstrap_dependencies.sh --prefix /opt/cbd-tools"

ENV PATH="/opt/cbd-tools/bin:/opt/conda/envs/cbdating/bin:${PATH}"

ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["source /opt/conda/etc/profile.d/conda.sh && conda activate cbdating && python scripts/preflight_check.py --methods all"]
