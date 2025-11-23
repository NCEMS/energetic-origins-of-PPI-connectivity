# Base: JupyterLab + mamba + notebook helpers
FROM jupyter/minimal-notebook:python-3.11

LABEL org.opencontainers.image.title="Annotated Network of the Yeast Interactome (ANYI)" \
      org.opencontainers.image.description="JupyterLab + annotated yeast PPI network"


# Use root for installs/permissions, then drop back to jovyan
USER root

# ---- Dependencies (conda/mamba) ----
# NOTE: Build context will be repo root; see build command below.
COPY --chown=${NB_UID}:${NB_GID} docker/environment.yml /tmp/environment.yml
RUN mamba env update -n base -f /tmp/environment.yml && \
    mamba clean -a -y && \
    fix-permissions "${CONDA_DIR}" && fix-permissions "/home/${NB_USER}"

# ---- Data (baked in) ----
# Your file lives at ../flatten/... relative to docker/, but since the build
# context is repo root, we copy from flatten/...
COPY --chown=${NB_UID}:${NB_GID} flatten/processed-data/20251023-s288c-annotated-PPI-net.pkl \
    /home/${NB_USER}/data/20251023-s288c-annotated-PPI-net.pkl

# Expose DB path as an environment variable
ENV PPI_DB_PATH=/home/${NB_USER}/data/20251023-s288c-annotated-PPI-net.pkl

# ---- Example notebooks (baked in to /examples) ----
COPY --chown=${NB_UID}:${NB_GID} docker/examples/ /home/${NB_USER}/examples/

# ---- Entry point (copies examples into /work on first run) ----
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

RUN mkdir -p /home/${NB_USER}/.local/share/jupyter/runtime && \
    chown -R ${NB_UID}:${NB_GID} /home/${NB_USER}/.local && \
    fix-permissions "/home/${NB_USER}"
ENV JUPYTER_RUNTIME_DIR=/tmp/jupyter

# Optional: named kernel in Jupyter
RUN python -m ipykernel install --user --name=ppi --display-name="Yeast PPI"

# JupyterLab by default; non-root user for runtime
ENV JUPYTER_ENABLE_LAB=yes
USER ${NB_UID}
EXPOSE 8888

# Start script copies examples → /work (if empty), then launches notebook server
CMD ["entrypoint.sh"]
