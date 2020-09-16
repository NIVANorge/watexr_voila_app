ARG BASE_CONTAINER=jupyter/base-notebook:6d42503c684f
FROM $BASE_CONTAINER

LABEL maintainer="James Sample <james.sample@niva.no>"

USER root

COPY requirements.txt /tmp/
RUN conda install -c conda-forge --quiet --yes --file /tmp/requirements.txt && \
    conda clean --all -f -y && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER && \
    rm -rf /tmp/*

COPY app /app/

USER $NB_UID

CMD ["voila", "/app/voila_app.ipynb", "--port", "8866", "--no-browser"] 