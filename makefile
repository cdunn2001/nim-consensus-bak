A=DW_banded.c   common.h      falcon.c      kmer_lookup.c poo.c

go: DW_banded.nim falcon.nim kmer_lookup.nim poo.nim common.nim

%.nim: %.c
	c2nim -o:$@ $<
%.nim: %.h
	c2nim -o:$@ $<
