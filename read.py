from bamread.src.bamread import bamread


def read_bam(f, mapq=0, required_flag=0, filter_flag=1540):

    return bamread(f, mapq, required_flag, filter_flag)
