# Snakemake notes

I had some decisions and notes I gathered along the way

## Pulling sequencing data using Snakemake

I wrote my own manual puller. However there is an option to use directly `snakemake`. Maybe I should use the native functionality. The problem was that our cluster had an older version without the functionality when I started this project, but it's updated now.

I want to pull data from NCBI and `snakemake` seems to have a module exactly for this. There is a function `snakemake.remote.NCBI` for pulling data from NCBI (here is its [documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#genbank-ncbi-entrez)). One would probably like to specify `--keep-remote` if this option will be used.

### other flags to consider

- `--cluster-status 'bjobs'` : allow snakemake to look at the status of jobs; ~this is not working on version of snakemake on Vital-it~ it was updated, but for some reason cluster status is not allowed even in version 11
- `--jobscript cluster_wrapper.sh`

## Snakemake tips

apparently I can produce a graph of the workflow (it's really pretty) :

```
snakemake --forceall --dag | dot -Tpng > dag1.png
```

Snakemake version of `make clean` (removed all downloaded and coputed data) is

```
rm $(snakemake --summary | tail -n+2 | cut -f1)
```

Snakemake version of GNU make dry run (show commands, don't execute them) is

```
snakemake -p --quiet -n
```
