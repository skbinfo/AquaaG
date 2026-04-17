#!/bin/bash
set -e

echo "======================================================="
echo "   AquaG Pipeline: Local Build & Initial Setup"
echo "======================================================="

# --- 1. Check for Docker ---
if ! command -v docker &> /dev/null; then
    echo "❌ ERROR: Docker is not installed or not in PATH."
    echo "Please install Docker natively on this machine (e.g., sudo apt install docker.io)"
    echo "and ensure your user is added to the 'docker' group: sudo usermod -aG docker \$USER"
    exit 1
fi
echo "✅ Docker detected: $(docker --version)"

# --- 2. Build the Main Pipeline Image Locally ---
echo "-------------------------------------------------------"
echo " Building Main Pipeline Image..."
echo "-------------------------------------------------------"
if [ ! -f "Dockerfile" ]; then
    echo "❌ ERROR: Dockerfile not found in the current directory."
    echo "Please run this script from the root of the cloned AquaG repository."
    exit 1
fi

docker build -t aquaag-pipeline:latest .
echo "✅ Main pipeline image built successfully."

# --- 3. GeneMark License Key (For Eukaryotes) ---
echo "-------------------------------------------------------"
echo " ACTION REQUIRED: GeneMark License Key (gm_key)"
echo " Get gm_key from: https://exon.gatech.edu/GeneMark/license_download.cgi"
echo "-------------------------------------------------------"

if [ -f "gm_key" ]; then
    echo "✅ 'gm_key' already found in current directory."
else
    while true; do
        read -p "Enter FULL path to your downloaded 'gm_key' file (or type 'skip'): " GM_KEY_PATH
        if [[ "$GM_KEY_PATH" == "skip" ]]; then
            echo "⏭️ Skipping gm_key setup. (You will need to place it in this folder before running Eukaryotic mode)."
            break
        elif [ -f "$GM_KEY_PATH" ]; then
            cp "$GM_KEY_PATH" ./gm_key
            echo "✅ License key copied to current directory as ./gm_key"
            break
        else
            echo "❌ ERROR: File not found."
        fi
    done
fi

# --- 4. EggNOG Database (Functional Annotation) ---
echo "-------------------------------------------------------"
echo " ACTION REQUIRED: EggNOG Database Setup"
echo "-------------------------------------------------------"
if [ -d "eggnog_db" ] && [ "$(ls -A eggnog_db 2>/dev/null)" ]; then
    echo "✅ EggNOG database already exists in ./eggnog_db"
else
    echo "Do you want to pre-download the 50GB EggNOG database now? (Required to use --run-func)"
    read -p "Enter [y/N]: " DOWNLOAD_EGGNOG

    if [[ "$DOWNLOAD_EGGNOG" =~ ^[Yy]$ ]]; then
        echo "-> Creating eggnog_db directory in $PWD..."
        mkdir -p "$PWD/eggnog_db"
        
        echo "-> Downloading databases via Docker (This will take a while)..."
        docker run --rm -it \
          -v "$PWD/eggnog_db":/eggnog_db \
          --entrypoint /bin/bash \
          quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_2 \
          -c "sed -i 's/eggnogdb.embl.de/eggnog5.embl.de/g' \$(which download_eggnog_data.py) && download_eggnog_data.py -P -y --data_dir /eggnog_db && chown -R $(id -u):$(id -g) /eggnog_db"
        
        echo "✅ EggNOG Database Setup Complete!"
    else
        echo "⏭️ Skipping EggNOG database download."
        echo "   (The pipeline will automatically download it if you use --run-func later)."
    fi
fi

echo "======================================================="
echo " You are ready to run the aquaag-pipeline command."
echo "======================================================="
