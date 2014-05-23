# ID Tidy

Daniel Standage, daniel.standage@gmail.com

The ``idtidy`` python module and associated scripts were written to facilitate minting clean, concise, and consistent IDs for sequences and annotations within the *Polistes dominula* genome project.

```bash
egrep -v '(expressed_sequence|protein|translated_nucleotide)_match' maker.gff3 | \
    gt gff3 -retainids -sort -tidy 2> gt.log | \
    python annot-ids.py --idfmt='Pdom%s-%05lu' -n --rnamap=rnaids.txt --dbxref=MAKER - | \
    python fix-regions.py scaffolds.fa > final-annot.gff3

python seq-ids.py rnaids.txt < maker-transcripts.fasta > final-transcripts.fasta
python seq-ids.py rnaids.txt < maker-proteins.fasta > final-proteins.fasta
```
