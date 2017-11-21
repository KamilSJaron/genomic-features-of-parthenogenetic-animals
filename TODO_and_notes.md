# TODO

map reads to corresponding sequences.

species by species
each sample with data to each reference available

for SP in YXXX;
    for REFERENCE in "$SP"[123]/genome.fasta.gz ;
        for SAMPLE in "$SP"[123]/reads_1.fq.gz ;
            map $REFERENCE $SAMPLE
        done
    done
done


# Notes

A springtail Holacanthella duospinosa is sequenced. Is it the closed sequenced relative of the asexual springtail genome?


