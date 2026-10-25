python python-scripts/make-Rosetta-commands.py generate \
    --nodes ../5-dG-calculations/processed-data/missing-s288c-step5.pkl \
            ../5-dG-calculations/processed-data/y2h-test-s288c-step5.pkl \
    --source-root .. \
    --bundle rosetta-transfer \
    --remote-bundle /path/on/supercomputer/rosetta-bundle \
    --relax-executable /path/on/supercomputer/relax.linuxgccrelease \
    --nstruct 10
