NIMFLAGS=-d:debug
NIMFLAGS=-d:release
NIMFLAGS+=--verbosity:2

run-main:
convert: DW_banded.nim falcon.nim kmer_lookup.nim poo.nim common.nim
run-%: %.exe
	./$*.exe
%.exe: %.nim
	nim ${NIMFLAGS} --out:$*.exe c $<
%.nim: %.c
	c2nim -o:$@ $<
%.nim: %.h
	c2nim -o:$@ $<
