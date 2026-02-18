FROM jupyter/minimal-notebook:python-3.11

LABEL org.opencontainers.image.title="ANnotated Yeast Interactome (ANYI)" \
      org.opencontainers.image.description="JupyterLab + annotated yeast PPI network browser tool"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      curl \
      ca-certificates \
      tar \
      gzip \
      unzip \
      bzip2 \
      xz-utils \
      findutils \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --chown=${NB_UID}:${NB_GID} docker/environment.yaml /tmp/environment.yaml
RUN mamba env update -n base -f /tmp/environment.yaml && \
    mamba clean -a -y && \
    fix-permissions "${CONDA_DIR}" && fix-permissions "/home/${NB_USER}"

COPY --chown=${NB_UID}:${NB_GID} docker/ANYI-browser/ /home/${NB_USER}/ANYI-browser/
ENV ANYI_BROWSER_NODES=/home/${NB_USER}/examples/ANYI-browser/data/nodes.pkl
ENV ANYI_BROWSER_EDGES=/home/${NB_USER}/examples/ANYI-browser/data/edges.csv

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

RUN mkdir -p /home/${NB_USER}/.local/share/jupyter/runtime && \
    mkdir -p /tmp/jupyter && \
    chown -R ${NB_UID}:${NB_GID} /home/${NB_USER}/.local /tmp/jupyter && \
    fix-permissions "/home/${NB_USER}"
ENV JUPYTER_RUNTIME_DIR=/tmp/jupyter

RUN python -m ipykernel install --user --name=ppi --display-name="Yeast PPI"

ENV JUPYTER_ENABLE_LAB=yes
USER ${NB_UID}
EXPOSE 8888
CMD ["entrypoint.sh"]
