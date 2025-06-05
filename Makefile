IMAGE_NAME = snakemake-msi
WORKDIR = /app

build:
	docker build -t $(IMAGE_NAME) .

run:
	docker run --rm -v "$$(pwd)":$(WORKDIR) -w $(WORKDIR) $(IMAGE_NAME)

run-debug:
	docker run --rm -v "$$(pwd)":$(WORKDIR) -w $(WORKDIR) $(IMAGE_NAME) \
	micromamba run -n snakemake-msi-env snakemake --cores 2 --rerun-incomplete

clean:
	rm -rf raw_data/*.txt

docker-clean:
	docker system prune -f --volumes
