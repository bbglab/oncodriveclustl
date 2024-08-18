FROM python:3.12-slim

# hadolint ignore=DL3042
RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,target=/project,rw \
    cd /project && pip install .

RUN oncodriveclustl --help

ENTRYPOINT [ "/usr/local/bin/oncodriveclustl" ]
