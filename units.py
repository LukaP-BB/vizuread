import matplotlib.pyplot as plt
from vizuread import plot_region, get_reads_from, parse_position, READ_SPACING

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(27,10))
f = "tests/child.bam"
plot_region(f, "chr22:23480074-23483000", ax=ax, piling=None)
plt.savefig("test.png")
plt.show()

# position = parse_position("chr22:1-100000000")
# for r in get_reads_from(f, *position) :
#     if "3D" in r.cigar :
#         print(r)
# c1 = "chr11"

# f = "T30989_realigned.fixed.recal.bam"
# c1 = "chr11"
# c2 = "chr14"
# s = "69,638,162"
# e = "69,639,433"

# reads = get_reads_from(f, c1, s, e, flags="-F 2")
# # s = "0"
# # e = "69,639,433"
# # for r in get_reads_from(f, c1, s, e) :
# #     if "3D" in r.cigar :
# #         print(r)

# # reads = get_reads_from(f, *position, samtools_command="samtools")
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12,5))

# for i, r in enumerate(reads) :
#     # print(r)
#     r.plot(ax, i)
# plt.show()

# print(is_forward(10))


# cigar = "35H8M1D32M"
# parse_cigar(cigar)