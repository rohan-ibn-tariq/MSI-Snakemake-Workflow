# FROM mambaorg/micromamba:1.4.5

# WORKDIR /app

# COPY environment/environment.yaml /tmp/environment.yaml

# RUN micromamba create -n snakemake-msi-env -f /tmp/environment.yaml -y && \
#     micromamba clean --all --yes

# COPY . /app/

# RUN chmod +x workflow/scripts/*.py || true
# RUN chmod +x workflow/scripts/*.awk || true

# SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

# ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]

# CMD ["micromamba", "run", "-n", "snakemake-msi-env", "snakemake", "--help"]


FROM mambaorg/micromamba:1.4.5

WORKDIR /app

# Copy environment file and create conda environment
COPY environment/environment.yaml /tmp/environment.yaml
RUN micromamba create -n snakemake-msi-env -f /tmp/environment.yaml -y && \
    micromamba clean --all --yes

# Copy all workflow files
COPY . /app/

# Make scripts executable
RUN chmod +x workflow/scripts/*.py || true
RUN chmod +x workflow/scripts/*.awk || true

# NO SHELL directive
# NO ENTRYPOINT directive

# Try to run snakemake directly
CMD ["snakemake", "--use-conda", "--cores", "4"]
