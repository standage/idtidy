#!/usr/bin/env python
import getopt, re, sys


# Supporting classes
class Entry():
  """
  Very simplistic representation of a GFF3 entry: a single line of a GFF3 file
  which may or may not be a complete representation of a genomic feature.
  """
  def __init__(self, line):
    self.fields = line.rstrip().split("\t")
    if not self.is_feature():
      self.line = line
      return
  
    self.seqid  = self.fields[0]
    self.source = self.fields[1]
    self.ftype  = self.fields[2]
    self.start  = self.fields[3]
    self.end    = self.fields[4]
    self.score  = self.fields[5]
    self.strand = self.fields[6]
    self.phase  = self.fields[7]
    self.attrs  = {}
    for keyvaluepair in self.fields[8].split(";"):
      key, value = keyvaluepair.split("=")
      if key in self.attrs:
        self.attrs[key] += ","+ value
      else:
        self.attrs[key] = value
  
  def is_feature(self):
    return len(self.fields) == 9

class GFF3IDMinter():
  """
  Class for minting IDs for gene and RNA features using the specified format.
  """
  def __init__(self, fp, idfmt = "%s%d"):
    """
    GFF3IDMinter constructor
    """
    self.oldids = {}
    self.scan_ids(fp)
    self.newids = {}
    self.genecount = 0
    self.mint_new_ids(idfmt)

  def scan_ids(self, instream):
    """
    Scan the GFF3-formatted data from the given input stream and store the
    mapping of gene IDs to RNA IDs
    """
    for line in instream:
      entry = Entry(line)
      if not entry.is_feature() or entry.ftype not in ("mRNA", "tRNA", "rRNA"):
        continue
      
      assert "ID" in entry.attrs, "Error: RNA features must have ID attributes"
      rnaid = entry.attrs["ID"]
      assert "Parent" in entry.attrs, ("Error: RNA features must be connected "
                                       "to a gene feature via its Parent "
                                       "attribute")
      geneid = entry.attrs["Parent"]
      
      if geneid not in self.oldids:
        self.oldids[geneid] = []
      self.oldids[geneid].append(entry)

  def mint_new_ids(self, idfmt):
    """
    Create new IDs for each gene and RNA feature.
    """
    for geneid, rnalist in self.oldids.iteritems():
      self.genecount += 1
      rnacount = 0
      newgeneid = idfmt % ("GENE", self.genecount)
      self.newids[geneid] = newgeneid
  
      for rnaentry in rnalist:
        rnacount += 1
        rnaid = rnaentry.attrs["ID"]
        rnatype = rnaentry.ftype.upper()
        rnafmt = idfmt+".%d"
        newrnaid = rnafmt % (rnatype, self.genecount, rnacount)
        self.newids[rnaid] = newrnaid
  
  def fix_line(self, line, dbxref = None):
    """
    If the given line from a GFF3 file is a feature, replace any gene or RNA IDs
    located therein.
    """
    entry = Entry(line)
    if not entry.is_feature():
      return line

    replaceid = "ID" in entry.attrs and entry.attrs["ID"] in self.newids
    replaceparent = "Parent" in entry.attrs and entry.attrs["Parent"] in self.newids 
    if not replaceid and not replaceparent:
      return line

    if replaceid:
      newid = self.newids[entry.attrs["ID"]]
      line = re.sub("ID=[^;]+", "ID="+newid, line)
      if dbxref:
        line += ";Dbxref=%s:%s" % (dbxref, entry.attrs["ID"])

    if replaceparent:
      newparent = self.newids[entry.attrs["Parent"]]
      line = re.sub("Parent=[^;]+", "Parent="+newparent, line)

    return line

  def write_genemap(self, fp):
    if not fp:
      return
    for geneid in self.oldids.keys():
      print >> fp, "%s\t%s" % (self.newids[geneid], geneid)

  def write_rnamap(self, fp):
    if not fp:
      return
    for rnalist in self.oldids.values():
      for entry in rnalist:
        rnaid = entry.attrs["ID"]
        print >> fp, "%s\t%s" % (self.newids[rnaid], rnaid)


# Other GFF3 processing functions
def strip_name(line):
  return re.sub(";*Name=[^;]+", "", line)

def strip_exon_id(line):
  if not re.search("\texon\t", line):
    return line
  return re.sub("ID=[^;]+;*", "", line)

def fix_cds_utr_id(line):
  labels = { "five_prime_UTR": "5putr", "three_prime_UTR": "3putr",
             "UTR": "utr", "CDS": "cds" }
  typematch = re.search("\t([^\t]*UTR|CDS)\t", line)
  if not typematch:
    return line

  ftype = typematch.group(1) 
  assert ftype in labels, "Error: unknown UTR type "+ ftype
  label = labels[typematch.group(1)]
  parentmatch = re.search("Parent=([^;]+)", line)
  assert parentmatch, "Error: CDS or UTR entry must have a Parent attribute"
  newid = "%s.%s" % (parentmatch.group(1), label)
  return re.sub("ID=[^;]+", "ID="+newid, line)


# Command line interface functions
def print_usage(outstream):
  usage = ("\n"
    "Create new IDs for gene and RNA features using the specified format, tested\n"
    "primarily on output from the Maker genome annotation pipeline.\n\n"
    "Usage: python annot-ids.py [options] annot.gff3\n"
    "  Options:\n"
    "    -f|--idfmt: STRING     printf-style format to use for creating new IDs;\n"
    "                           must accept a string (%s) and long unsigned\n"
    "                           integer (%d or %lu); default is '%s%d'\n"
    "    -g|--genemap: FILE     write the correspondence of new gene IDs to old\n"
    "                           gene IDs to the given file\n"
    "    -h|--help              print this help message and exit\n"
    "    -n|--stripnames        remove names from all features; particularly useful\n"
    "                           when name attributes are uninformative, as is the\n"
    "                           case with Maker\n"
    "    -o|--outfile: FILE     file to which output will be written; default is\n"
    "                           the terminal (standard output)\n"
    "    -r|--rnamap: FILE      write the correspondence of new RNA IDs to old\n"
    "                           RNA IDs to the given file\n"
    "    -x|--dbxref: STRING    use the 'Dbxref' attribute (database cross\n"
    "                           reference) to store a copy of any ID that is\n"
    "                           replaced; the argument provided will serve as the\n"
    "                           Dbxref key\n")
  print >> outstream, usage

def parse_options(argv):
  params = {
    "dbxref": None,
    "genemap": None,
    "idfmt": "%s%d",
    "outfile": sys.stdout,
    "rnamap": None,
    "stripnames": False,
    "infile": None
  }
  optstr = "f:g:hno:r:x:"
  longopts = ["idfmt=", "genemap=", "help", "stripnames", "outfile=", "rnamap=",
              "dbxref="]
  (options, args) = getopt.getopt(argv[1:], optstr, longopts)
  for key, value in options:
    if key in ("-f", "--idfmt"):
      params["idfmt"] = value
    elif key in ("-g", "--genemap"):
      params["genemap"] = open(value, "w")
    elif key in ("-h", "--help"):
      print_usage(sys.stdout)
      sys.exit(0)
    elif key in ("-n", "--stripnames"):
      params["stripnames"] = True
    elif key in ("-o", "--outfile"):
      params["outfile"] = open(value, "w")
    elif key in ("-r", "--rnamap"):
      params["rnamap"] = open(value, "w")
    elif key in ("-x", "--dbxref"):
      params["dbxref"] = value
    else:
      print_usage(sys.stderr)
      assert False, "unsupported option '%s'" % key

  if len(args) != 1:
    print_usage(sys.stderr)
    assert False, "expected 1 argument, %d provided" % len(args)
  if args[0] == "-":
    params["infile"] = sys.stdin
  else:
    params["infile"] = open(args[0])

  return params

if __name__ == "__main__":
  params = parse_options(sys.argv)

  # Load data into memory; increases runtime and memory consumption, but
  # required for correct handling of stdin, named pipes, or process substitutions
  indata = []
  for line in params["infile"]:
    indata.append(line)

  minter = GFF3IDMinter(indata, params["idfmt"])
  for line in indata:
    line = line.rstrip()
    line = minter.fix_line(line, params["dbxref"])
    if params["stripnames"]:
      line = strip_name(line)
    line = strip_exon_id(line)
    line = fix_cds_utr_id(line)
    print >> params["outfile"], line

  minter.write_genemap(params["genemap"])
  minter.write_rnamap(params["rnamap"])
