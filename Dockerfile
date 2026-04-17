FROM condaforge/miniforge3:latest

# Install Docker CLI so this container can talk to the host's Docker daemon
# (Required to spin up BRAKER3 and EggNOG-mapper containers)
RUN apt-get update && \
    apt-get install -y docker.io && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

# Copy the main conda environment file into the container FIRST
# Note: Ensure eggnog-mapper is REMOVED from your environment_main.yml
COPY environment_main.yml .

# 1. Limit concurrent downloads to prevent bandwidth starvation (fetch_threads 2)
# 2. Try to create the env. If the network fails, use 'env update' in a loop to resume!
# 3. Clean up the package cache to keep the final image size down
RUN conda config --set fetch_threads 2 && \
    conda config --set remote_read_timeout_secs 1200 && \
    conda config --set remote_connect_timeout_secs 100 && \
    conda config --set remote_max_retries 10 && \
    mamba env create -f environment_main.yml || \
    for i in 1 2 3 4 5 6 7 8 9 10; do \
        echo "Network timeout hit. Resuming downloads (Attempt $i/10)..." && \
        mamba env update -n assembly_tool_env_new -f environment_main.yml && break || sleep 5; \
    done && \
    mamba clean -afy

# Copy your pipeline script and the dynamic configuration script
COPY AquaaG.py .
COPY generate_dashboard.py .
COPY entrypoint.sh .

# Make both scripts executable
RUN chmod +x AquaaG.py entrypoint.sh

# Create a default directory for Augustus config to prevent path errors in BRAKER
RUN mkdir -p /augustus_config

# Set the entrypoint to the bash wrapper script
ENTRYPOINT ["/app/entrypoint.sh"]
