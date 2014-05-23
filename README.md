# ID Tidy
Daniel Standage, <daniel.standage@gmail.com>

The ``idtidy`` python module and associated scripts were written to facilitate minting clean, concise, and consistent IDs for sequences and annotations within the *Polistes dominula* genome project.
A complete description of settings for the ``annot-ids.py`` program are available by executing ``python annot-ids.py --help`` on the command line.
The example below demonstrates the use of the other two scripts in context.

```bash
python annot-ids.py --idfmt='Pdom%sr1.0-%05lu' \
                    --genemap=newgeneids.txt \
                    --rnamap=newrnaids.txt \
                    --dbxref=MAKER \
                    --outfile=annot-newids.gff3
                    annot.gff3
python fix-regions.py scaffolds.fa < annot-newids.gff3 > annot-fixed.gff3
python seqids.py newrnaids.txt < annot-transcripts.fasta > annot-trnascripts-newids.fasta
```


