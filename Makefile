# Basic variables
IMAGE_NAME = annotated-network-yeast-interactome
DOCKERFILE = docker/Dockerfile
TOKEN ?= letmein
NOTEBOOK_PORT ?= 8888

# Default target
.PHONY: help
help:
	@echo "Usage:"
	@echo "  make build            - Build the Docker image"
	@echo "  make run              - Run container (ephemeral, in-memory)"
	@echo "  make run-persistent   - Run with ./user_notebooks mounted for persistence"
	@echo "  make clean            - Remove stopped containers and dangling images"

# Build the image from repo root
.PHONY: build
build:
	docker build -t $(IMAGE_NAME) -f $(DOCKERFILE) .

# Run without mounting (changes lost on exit)
.PHONY: run
run:
	docker run --rm -p $(NOTEBOOK_PORT):8888 \
		-e JUPYTER_TOKEN=$(TOKEN) \
		$(IMAGE_NAME)

# Run with persistent host directory for user notebooks
.PHONY: run-persistent
run-persistent:
	mkdir -p user_notebooks
	docker run --rm -p $(NOTEBOOK_PORT):8888 \
		-v $$PWD/user_notebooks:/home/jovyan/work \
		-e JUPYTER_TOKEN=$(TOKEN) \
		$(IMAGE_NAME)

# Clean dangling images and containers
.PHONY: clean
clean:
	docker system prune -f
