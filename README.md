AGL
===

Assign Germ Line
----------------

A germline assignment program for antibody protein sequences

Installation
------------

First you must install `BiopLib` (https://github.com/ACRMGroup/bioplib/).

Then for most people, simply type:
```
./install.pl
```
to install in `~/bin`, or type:
```
./install.pl dirname
```
to install in `dirname`.

This will compile and install the code and place data files in a
sub-directory of where the executable is installed.

You can also move the data files elsewhere and set the environment
variable, `AGLDATADIR` to point to the directory containing the data
files (or use `-d` on the `agl` command line to specify the location
of the files)..

Philosophical problems
----------------------

There is a philosophical problem with germline assignment.

Take this *Homo sapiens* VH sequence from INN 11258

```
EVQLVESGPGLVKPSETLSLTCTVSGFSLTRFGVHWVRQPPGKGLEWLGV
IWRGGSTDYNAAFVSRLTISKDNSKNQVSLKLSSVTAADTAVYYCSNHYY
GSSDYALDNWGQGTLVTVSS
```

- It matches IGHV4-4\*08 with 70/97 residues (72.16%) and
- it matches IGHV4-59\*03 with 70/96 residues (72.92%)

The second one has a higher percent ID, but the first has better
coverage (the germline sequence covers more residues - the R is
missing from the end of the second one's sequence).

IMGT-domainGapAlign selects IGHV4-4\*08

Bizarrely, however, domainGapAlign shows them both as 72.9% sequence
and says there is an overlap of 96 residues for both (this is not the
case)! 

It also lists IGHV4-59\*01, IGHV4-59\*08, IGHV4-59\*13, IGHV4-59\*02
as having the same percent identity (72.9%) and overlap (96 residues)
and all but IGHV4-59\*02 have the same Smith Waterman score (462;
IGHV4-59\*02 scores 461). I assume that this is calculated from the
DNA information - i.e. how many base changes are needed for an amino
acid change - but this is irrelevant for antibody drugs as these will
have been codon optimized to improve expression.

Actual Known Problems
---------------------

- 11351L gives VL IGKV3-11\*01 instead of IGKV3D-11\*02 (longer but a lot worse)

- 11366H gives JH IGHJ3\*01 instead of IGHJ1\*01

- 11390L gives VL IGLV2-14\*01 instead of IGLV2-14\*03 (which is longer and better)

```
QTVVTQEPSLTVSPGGTVTLTCGSSTGAVTSGNYPNWVQQKPGQAPRGLI
GGTKFLAPGTPARFSGSLLGGKAALTLSGVQPEDEAEYYCVLWYSNRWVF
GGGTKLTVL
```
- gives JL IGLJ3\*01 instead of IGLJ2\*01 (which is identical)

- 12094L gives JL IGLJ3\*01 instead of IGLJ2\*01 (which is identical)

- 11881H gives JH IGHJ4\*01 instead of IGHJ1\*01

- 12124H gives JH IGHJ4\*01 instead of IGHJ1\*01

- 11672L gives VL IGLV2-14\*02 instead of IGLV2-23\*02 (much longer!)

- `BUGS/fails.faa` is a heavy chain that fails altogether (it seems to
  have a deletion in HFR3 and a low seqid, but if we tell `agl` this
  is human heavy chain it should still report a result)

- 11878L doesn't find a JL at all (should be *Homo sapiens* IGKJ5\*01)

- 12056H gives IGHJ4*01 rather than IGHJ6*01 though the latter is
  longer with the same number of mismatches

- 11938H (trastuzumab) gives IGHV3-66*01 rather than IGHV3-23*04 The
  latter scores better but has an indel. On the other hand 12079H does
  give IGHV3023*04 (but BLAST gives IGHV3-66*01)

- 12050H gives CH1,CH2 and CH1\*01 and IGHG1\*01 as IGHG1\*07 - they should
  all be IGHG1\*07

- 12011L doesn't assign J - should be IGKJ5\*01 (note has deletions)

- 11943H gives JH  IGHJ4\*01 instead of IGHJ5\*01 (longer)

- 11857L gives the Distal sequence IGKV1D-13\*01 rather than the normal
  one since there is no \*01 for IGKV1-13 (i.e. \*01 is taking
  precedence over not using D)

- 12145H gives JH IGHJ4\*01 instead of IGHJ1\*01 1\*01 is the same length
  and has the same score so should be selected on the basis of lower number

- 21319H[5] gives VH IGHV3-74\*01 instead of IGHV3-66\*01 (which has an insert
  in a CDR). Similarly trastuzumab gives the wrong VH because of an insertion

- 12528H gives JH IGHJ4\*01 instead of IGHJ2\*01 which is the same length and
  had better seqid IGHJ4\*01: 12/14 IGHJ2*01: 13/14
  This is because the translation for IGHJ2\*01 is 2 residues longer at the
  N-terminus
```  
  >J-REGION_IGHJ2*01_F1_Homo sapiens
  YWYFDLWGRGTLVTVSS

  >J-REGION_IGHJ4*01_F2_Homo sapiens
    YFDYWGQGTLVTVSS
```

- 12552H gives JH IGHJ4\*01 instead of IGHJ1\*01. Looking at the alignment
  they seem to have the same number of mismatches, but 4\*01 scores
  more highly. I'm guessing that it is actually longer. I think we need to
  trim tails from the sequences before calculating the score.

- 12545H gives mouse VH IGHV1-53*01 instead of IGHV1S81*02 which is one
  mutation better.

- 12421H gives JH IGHJ3\*01 instead of IGHJ6*01. The latter is longer
  and has the same number of mismatches.

Data and setup
--------------

The code uses germline DNA data from IMGT. The `data` directory of the
distribtion contains these Fasta files containing DNA sequences for
the germline gene segments obtained from
http://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgtrefseq.fasta

`makedb.pl` translates these to create files in `./share/agl/data`.

By default, AGL will expect this directory to exist under the
directory in which the agl executable is found (the `install.sh`
script will handle all this). Alternatively, the environment variable,
`AGLDATADIR` can point this data directory (or you can use `-d` on the
`agl` command line).

**Note that commercial use of the DNA data requires a licence from IMGT.**
