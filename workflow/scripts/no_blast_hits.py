import skbio
import pandas as pd

df = skbio.io.read(
    snakemake.input[1],
    format="blast+7",
    into=pd.DataFrame,
)

with open(snakemake.output[0], "w+") as f:
    for seq in skbio.read(snakemake.input[0], "fasta"):
        if not len(
            df.loc[
                (df.qseqid == seq.metadata["id"])
                & (df.evalue < snakemake.params["max_evalue"])
            ]
        ):
            seq.write(f)
