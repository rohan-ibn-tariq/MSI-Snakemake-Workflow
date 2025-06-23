FROM mambaorg/micromamba:1.4.5

WORKDIR /app

COPY . /app/

COPY environment/environment.yaml /app/environment.yaml

RUN micromamba create -n snakemake-msi-env -f environment/environment.yaml -y && \
    micromamba clean --all --yes

RUN micromamba install -n snakemake-msi-env mamba -y

SHELL ["/bin/bash", "-c"]

CMD micromamba run -n snakemake-msi-env \
    snakemake --use-conda --cores 4
