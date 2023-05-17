from pathlib import Path
from typing import Tuple

import pandas as pd

from bamread.src.bamread import _bamread, _bamread_all  # type: ignore


def read_bam(
    f: Path, mapq: int = 0, required_flag: int = 0, filter_flag: int = 1540
) -> pd.DataFrame:
    chromosomes, starts, ends, strands, flags, chrmap = _bamread(
        f, mapq, required_flag, filter_flag
    )

    chromosomes, ends, flags, starts, strands = _create_series(
        chrmap, chromosomes, ends, flags, starts, strands
    )

    return pd.DataFrame(
        {
            "Chromosome": chromosomes,
            "Start": starts,
            "End": ends,
            "Strand": strands,
            "Flag": flags,
        }
    )


def _create_series(
    chrmap, chromosomes, ends, flags, starts, strands
) -> Tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    chromosomes = pd.Series(chromosomes).replace(chrmap).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).replace({16: "-", 0: "+"}).astype("category")
    flags = pd.Series(flags)
    return chromosomes, ends, flags, starts, strands


def read_bam_full(
    f: Path, mapq: int = 0, required_flag: int = 0, filter_flag: int = 1540
) -> pd.DataFrame:
    (
        chromosomes,
        starts,
        ends,
        strands,
        flags,
        chrmap,
        qstarts,
        qends,
        query_names,
        query_sequences,
        cigarstrings,
        query_qualities,
    ) = _bamread_all(f, mapq, required_flag, filter_flag)

    chromosomes, ends, flags, starts, strands = _create_series(
        chrmap, chromosomes, ends, flags, starts, strands
    )
    qstarts = pd.Series(qstarts)
    qends = pd.Series(qends)
    query_names = pd.Series(query_names)
    query_sequences = pd.Series(query_sequences)
    cigarstrings = pd.Series(cigarstrings)
    query_qualities = pd.Series(query_qualities)

    return pd.DataFrame(
        {
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
            "Quality": query_qualities,
        }
    )
