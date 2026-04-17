#!/bin/bash
set -e

# 1. Grab environment variables
EMAIL="${EMAIL:-researcher@example.com}"
TYPE="${TYPE:-PK}"
GROUP="${GROUP:-bacteria}"
ORGANISM="${ORGANISM:-}"
THREADS="${THREADS:-8}"
BUSCO="${BUSCO:-bacteria_odb10}"
PROKKA_KINGDOM="${PROKKA_KINGDOM:-Bacteria}"

# --- NEW: Auto-Detect GeneMark Key ---
# If the user dropped their gm_key in the current folder, automatically wire it up!
if [ -f "gm_key" ]; then
    export BRAKER_GM_KEY="$PWD/gm_key"
    echo "✅ Auto-detected gm_key in current directory!"
fi
# -------------------------------------

# 2. Set our temporary config path
CONFIG_PATH="/tmp/auto_config.yaml"

# 3. Write the common base configuration
cat <<EOF > $CONFIG_PATH
email: "$EMAIL"
organism_type: "$TYPE"
group: "$GROUP"
organism: "$ORGANISM"
quast_threads: $THREADS
busco_lineage: "$BUSCO"
busco_params:
  cpu: $THREADS
EOF

# 4. Append specific parameters based on EK or PK type
if [ "$TYPE" = "EK" ]; then
cat <<EOF >> $CONFIG_PATH
braker_params:
  cpus: $THREADS
EOF
else
cat <<EOF >> $CONFIG_PATH
prokka_kingdom: "$PROKKA_KINGDOM"
prokka_params:
  cpus: $THREADS
EOF
fi

# Print it out so the user can verify the settings
echo "=== Auto-generated Config ==="
cat $CONFIG_PATH
echo "============================="

# 5. Execute the script pointing to the /tmp/ config file
exec conda run --no-capture-output -n assembly_tool_env_new python /app/AquaaG.py -c $CONFIG_PATH "$@"
