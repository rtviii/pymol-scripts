reinitialize

fetch 1F9J, async=0
fetch 1YX5, async=# start a python block
python

# get two structures
cmd.fetch('2xwu 2x19', async=0)

# align and get raw alignment
cmd.align('/2xwu//B//CA', '/2x19//B//CA', cycles=0, transform=0, object='aln')
raw_aln = cmd.get_raw_alignment('aln')

# print residue pairs (atom index)
for idx1, idx2 in raw_aln:
    print('%s`%d -> %s`%d' % tuple(idx1 + idx2))

#end the python block
python end