IMAGE_NAME = snakemake-msi
CONTAINER_WORKDIR = /app
HOST_WORKDIR = $(shell pwd)

.PHONY: build run run-debug dry-run shell clean docker-clean test

# Build the Docker image
build:
	docker build -t $(IMAGE_NAME) .

# Run the full workflow
run: build
	@mkdir -p .snakemake
	docker run --rm \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/workflow/resources:$(CONTAINER_WORKDIR)/workflow/resources \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --use-conda --conda-frontend mamba --cores 4

# Debug run with rerun-incomplete
run-debug: build
	@mkdir -p .snakemake
	docker run --rm \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/workflow/resources:$(CONTAINER_WORKDIR)/workflow/resources \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --use-conda --conda-frontend mamba --cores 2 --rerun-incomplete

# Dry run to test workflow
dry-run: build
	@mkdir -p .snakemake
	docker run --rm \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --dry-run --cores 1

# Test with specific rule
test-rule: build
	@mkdir -p .snakemake
	docker run --rm \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --use-conda --conda-frontend mamba --cores 2 $(RULE)

# Interactive shell for debugging
shell: build
	docker run --rm -it \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/workflow/resources:$(CONTAINER_WORKDIR)/workflow/resources \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		bash

# Force rerun specific rule
force-rule: build
	@mkdir -p .snakemake
	docker run --rm \
		-v "$(HOST_WORKDIR)"/results:$(CONTAINER_WORKDIR)/results \
		-v "$(HOST_WORKDIR)"/workflow/logs:$(CONTAINER_WORKDIR)/workflow/logs \
		-v "$(HOST_WORKDIR)"/.snakemake:$(CONTAINER_WORKDIR)/.snakemake \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --use-conda --conda-frontend mamba --cores 2 --forcerun $(RULE)

# Clean up generated files
clean:
	rm -rf results/*
	rm -rf workflow/logs/*
	rm -rf .snakemake/

# Clean Docker artifacts
docker-clean:
	docker rmi $(IMAGE_NAME) || true
	docker system prune -f --volumes

# Show workflow DAG
dag: build
	docker run --rm \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --dag | dot -Tsvg > workflow_dag.svg

# List all targets
list-rules: build
	docker run --rm \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		snakemake --list

# Check environment
check-env: build
	docker run --rm \
		$(IMAGE_NAME) \
		micromamba run -n snakemake-msi-env \
		python --version && \
		snakemake --version
