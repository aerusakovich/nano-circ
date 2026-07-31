# Build from simulator/new/:  docker build -f containers/Dockerfile -t nanocirc-sim .
# Can also be converted straight to an Apptainer/Singularity image with:
#   apptainer build nanocirc-sim.sif docker-daemon://nanocirc-sim:latest
FROM condaforge/miniforge3:24.3.0-0

COPY environment.yml /opt/nanocirc-sim/environment.yml
RUN mamba env update -n base -f /opt/nanocirc-sim/environment.yml \
    && mamba clean -a -y

COPY pyproject.toml /opt/nanocirc-sim/pyproject.toml
COPY nanocirc_sim /opt/nanocirc-sim/nanocirc_sim
COPY configs /opt/nanocirc-sim/configs
RUN pip install --no-deps /opt/nanocirc-sim

ENTRYPOINT ["nanocirc-sim"]
CMD ["--help"]
