FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-numpy python3-scipy python3-requests python3-pip \
    quantum-espresso ca-certificates curl git time && \
    rm -rf /var/lib/apt/lists/*
RUN useradd -m -u 1000 agent && mkdir -p /app /logs/artifacts && \
    chown -R agent:agent /app /logs/artifacts
WORKDIR /app
CMD ["sleep", "infinity"]
