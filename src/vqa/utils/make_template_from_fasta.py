fasta_template = "../data/data/hcov_global_2023-11-16_09-28/template.fasta"

# read fasta file
with open(fasta_template, "r") as f:
    lines = f.readlines()

id = lines[0].strip()
seq = lines[1].strip()

# create template
# split sequence into 1000 bp chunks
chunks = [seq[i:i+1000] for i in range(0, len(seq), 1000)]

# template is in the format:
# >id:start_pos:end_pos
# seq
template = ""
for i, chunk in enumerate(chunks):
    template += f">{id}:{i*1000}:{i*1000 + len(chunk)}\n"
    template += chunk + "\n"


# write template to file
with open("../data/data/hcov_global_2023-11-16_09-28/template_1000bp.template", "w") as f:
    f.write(template)
