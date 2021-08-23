import pandas as pd

from bamread.src.bamread import _bamread, _bamread_all


def read_bam(f, mapq=0, required_flag=0, filter_flag=1540):

    chromosomes, starts, ends, strands, flags, chrmap = _bamread(
        f, mapq, required_flag, filter_flag)

    chromosomes = pd.Series(chromosomes).replace(chrmap).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).replace({16: "-", 0: "+"}).astype("category")
    flags = pd.Series(flags)

    return pd.DataFrame({
        "Chromosome": chromosomes,
        "Start": starts,
        "End": ends,
        "Strand": strands,
        "Flag": flags
    })


def read_bam_full(f, mapq=0, required_flag=0, filter_flag=1540):

    chromosomes, starts, ends, strands, flags, chrmap, qstarts, qends, query_names, query_sequences, cigarstrings, query_qualities = _bamread_all(
        f, mapq, required_flag, filter_flag)

    chromosomes = pd.Series(chromosomes).replace(chrmap).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).replace({16: "-", 0: "+"}).astype("category")
    flags = pd.Series(flags)
    qstarts = pd.Series(qstarts)
    qends = pd.Series(qends)
    query_names = pd.Series(query_names)
    query_sequences = pd.Series(query_sequences)
    cigarstrings = pd.Series(cigarstrings)
    query_qualities = pd.Series(query_qualities)

    return pd.DataFrame({
        "Chromosome": chromosomes,
        "Start": starts,
        "End": ends,
        "Strand": strands,
        "Flag": flags,
        "QueryStart": qstarts,
        "QueryEnd": qends,
        "QuerySequence": query_sequences,
        "Name": query_names,
        "Cigar": cigarstrings,
        "Quality": query_qualities
    })
