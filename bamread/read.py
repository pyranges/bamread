import pandas as pd

from bamread.src.bamread import _bamread


def read_bam(f, mapq=0, required_flag=0, filter_flag=1540):

    chromosomes, starts, ends, strands, flags, chrmap = _bamread(
        f, mapq, required_flag, filter_flag)

    chromosomes = pd.Series(chromosomes).replace(chrmap).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).replace({16: "+", 0: "-"}).astype("category")
    flags = pd.Series(flags)

    return pd.DataFrame({
        "Chromosome": chromosomes,
        "Start": starts,
        "End": ends,
        "Strand": strands,
        "Flag": flags
    })
