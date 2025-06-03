FROM mambaorg/micromamba:1.4.5

WORKDIR /app

COPY environment/environment.yaml /app/environment.yaml
COPY Snakefile /app/Snakefile
COPY raw_data /app/raw_data

RUN micromamba create -n snakemake-msi-env -f environment.yaml -y && \
    micromamba clean --all --yes

SHELL ["/bin/bash", "-c"]

CMD ["micromamba", "run", "-n", "snakemake-msi-env", "snakemake", "--cores", "2"]
