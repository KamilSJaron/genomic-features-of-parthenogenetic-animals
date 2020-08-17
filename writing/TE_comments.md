
00:13 15 Aug
this is puzzling - is this still the case? what data did you use? I've used dnaPipeTE on the Ar Illumina reads and get similar results as to Av: 1.4% and 1.8% for Ar and Av respectively, counting only interspersed TEs.. these are both much lower than the values in Fig 4, which is strange (let's chat about this!)
00:14 15 Aug
it might depend on what is your denominator on the X axis of Fig 4?
00:33 15 Aug
I also think there is an interesting aspect of dnaPipeTE when it comes to bdelloids anyway: as I understand it, dnaPipe works by identifying repeats as overrepresented regions in sequencing data, relative to the genomic (non-repeated) background, then assembles them (Trinity) and annotates them (RepeatMasker, BLAST etc). The problem with bdelloids is they contain a lot of low/single copy TEs, that might not be repeated above the tetraploidy level anyway, and I think things like this are possibly missed by dnaPipeTE. For bdelloids, a better approach (generally speaking, maybe not here for reasons below) might be to use RM on the assembled scaffolds and dnaPipe on the reads, and then the true value is approx the intersection of the two results. Using this approach, I get 4.21% for Ar and 5.76% for Av, which is about right for Av at least (according to the work of Irina Arkhipova, the maestro on bdelloid TEs). I'm not suggesting you redo the analysis (apart from checking that low Ar value) as it would be inconsistent across the taxa, but something to think about in the interpretation.
10:58 19 Aug
we used the reads with DNApipeTE as well. this is strange...
And yes, with the genome structure this might be a bit problematic.
I think we ran the Ar analyses a couple of times to double check. hm.
20:04 19 Aug
This comment will be unredeable soon. I suggest to move the discussion here:
https://github.com/KamilSJaron/genomic-features-of-asexual-animals/issues/22