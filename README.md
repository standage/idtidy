# Maker ID Mint

Daniel Standage, daniel.standage@gmail.com

The [Maker annotation pipeline](http://www.yandell-lab.org/software/maker.html) uses long, verbose IDs for the genomic features it annotates.
No doubt there is a purpose, at least within the Maker framework, for formatting the IDs the way they do.
However, when it comes to distribution these ID are pretty unappealing.

This is a small collection of scripts for replacing IDs for gene and RNA features in a GFF3 file.
These scripts were designed for use with any GFF3-formatted annotations, but they have only be tested on Maker-produced GFF3.

Here is a sample invocation of the scripts.

```bash
# The line containing 'gt gff3' can be removed if GenomeTools is not available on your system.
egrep -v '(expressed_sequence|protein|translated_nucleotide)_match' maker.gff3 | \
    gt gff3 -retainids -sort -tidy 2> gt.log | \
    python annot-ids.py --idfmt='Pdom%s-%05lu' -n --rnamap=rnaids.txt --dbxref=MAKER - | \
    python fix-regions.py scaffolds.fa > final-annot.gff3

python seq-ids.py rnaids.txt < maker-transcripts.fasta > final-transcripts.fasta
python seq-ids.py rnaids.txt < maker-proteins.fasta > final-proteins.fasta
```
