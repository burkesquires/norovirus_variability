from Bio import SeqIO

input_handle = open("noro-var.fasta", "rU")
output_handle = open("noro-var.faa", "w")

sequences = SeqIO.parse(input_handle, "fasta")
for sequence in sequences:
    sequence.seq.translate()
count = SeqIO.write(sequences, output_handle, "fasta")

output_handle.close()
input_handle.close()
print
"Coverted %i records" % count
