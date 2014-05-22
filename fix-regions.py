#!/usr/bin/env python
import re, sys

def parse_fasta(fp):
  """
  Stolen shamelessly from http://stackoverflow.com/a/7655072/459780.
  """
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, ''.join(seq))
      name, seq = line, []
    else:
      seq.append(line)
  if name: yield (name, ''.join(seq))

if __name__ == "__main__":
  usage = ("\nFix '##sequence-region' pragmas in the given GFF3 file.\n"
           "Usage: python fix-regions.py seqs.fasta < annot.gff3 > annot-fixed.gff3\n")

  if len(sys.argv) != 2:
    print >> sys.stderr, usage
  
  seqlengths = {}
  with open(sys.argv[1], 'r') as fp:
    for defline, seq in parse_fasta(fp):
      seqid = re.search("^>(\S+)", defline).group(1)
      seqlengths[seqid] = len(seq)

  for line in sys.stdin:
    line = line.rstrip()
    seqregmatch = re.search("^##sequence-region\s+(\S+)", line)
    if seqregmatch:
      seqid = seqregmatch.group(1)
      assert seqid in seqlengths, "sequence '%s' not found" % seqid
      line = "##sequence-region   %s 1 %lu" % (seqid, seqlengths[seqid])
    print line
  