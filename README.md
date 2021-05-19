AGL
===

Assign Germ Line
----------------

A germline assignment program for antibody protein sequences

Setup
-----

The `data` directory contains fasta files containing DNA sequences for
the germline gene segments.

`makedb.pl` translates these to create files in `mydata`.

The environment variable, `AGLDATADIR` must point to the `mydata`
directory (or use `-d` on the `agl` command line.


Philosophical problems
----------------------

There is a philosophical problem with germline assignment.

Take this Homo sapiens VH sequence from 11258

```
EVQLVESGPGLVKPSETLSLTCTVSGFSLTRFGVHWVRQPPGKGLEWLGV
IWRGGSTDYNAAFVSRLTISKDNSKNQVSLKLSSVTAADTAVYYCSNHYY
GSSDYALDNWGQGTLVTVSS
```

- It matches IGHV4-4*08 with 70/97 residues (72.16%) and
- it matches IGHV4-59*03 with 70/96 residues (72.92%)

The second one has a higher, yet the first has better coverage (the
germline sequence covers more residues - the R is missing from the
end of the second one's sequence).

domainGapAlign selects IGHV4-4*08

Bizarrely, however, domainGapAlign shows them both as 72.9% sequence
and says there is an overlap of 96 residues for both (this is not the
case)! 

It also lists IGHV4-59*01, IGHV4-59*08, IGHV4-59*13, IGHV4-59*02 as
having the same percent identity (72.9%) and overlap (96 residues) and
all but IGHV4-59*02 have the same Smith Waterman score (462;
IGHV4-59*02 scores 461). I assume that this is calculated from the DNA
information - i.e. how many base changes are needed for an amino acid
change. 

Actual Problems
---------------

- 11351L gives VL IGKV3-11*01 instead of IGKV3D-11*02 (longer but a lot worse)

- 11366H gives JH IGHJ3*01 instead of IGHJ1*01

- 11390L gives VL IGLV2-14*01 instead of IGLV2-14*03 (which is longer and better)

```
QTVVTQEPSLTVSPGGTVTLTCGSSTGAVTSGNYPNWVQQKPGQAPRGLI
GGTKFLAPGTPARFSGSLLGGKAALTLSGVQPEDEAEYYCVLWYSNRWVF
GGGTKLTVL
```
- gives JL IGLJ3*01 instead of IGLJ2*01 (which is identical)

- 12094L gives JL IGLJ3*01 instead of IGLJ2*01 (which is identical)

- 11881H gives JH IGHJ4*01 instead of IGHJ1*01

- 11672L gives VL IGLV2-14*02 instead of IGLV2-23*02 (much longer!)

- fails.faa is a heavy chain and fails altogether (seems to have a
  deletion in HFR3 and a low seqid, but if we tell agl this is human
  heavy chain it should still report a result)

- 11878L doesn't find a JL at all (should be Homo sapiens IGKJ5*01)
