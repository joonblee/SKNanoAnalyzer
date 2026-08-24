apptainer exec \
  /data6/Users/joonblee/SKNanoContainer/SKNANOAnalyzer_v13.sif \
  bash -lc '
    export PATH=/opt/conda/bin:$PATH
    export MAMBA_ROOT_PREFIX=/opt/conda

    eval "$(micromamba shell hook -s bash)"
    micromamba activate Nano

    cd /data6/Users/joonblee/SKNanoAnalyzer_v13
    source setup.sh

    ./scripts/build.sh
  '
