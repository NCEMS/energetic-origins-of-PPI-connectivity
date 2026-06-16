#snakemake -c all --use-conda --conda-frontend conda --snakefile 2-network-centrality/Snakefile --configfile config-files/union-athaliana.config
#snakemake -c all --use-conda --conda-frontend conda --snakefile 2-network-centrality/Snakefile --configfile config-files/in-vitro-athaliana.config
#snakemake -c all --use-conda --conda-frontend conda --snakefile 2-network-centrality/Snakefile --configfile config-files/in-vivo-athaliana.config
snakemake -c all --use-conda --conda-frontend conda --snakefile 2-network-centrality/Snakefile --configfile config-files/intersection-in-vitro-athaliana.config
snakemake -c all --use-conda --conda-frontend conda --snakefile 2-network-centrality/Snakefile --configfile config-files/intersection-in-vivo-athaliana.config
