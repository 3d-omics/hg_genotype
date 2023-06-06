def get_reads(wildcards):
    """Get reads for a sample and library."""
    forward_, reverse_ = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ][["forward", "reverse"]].values[0]
    return forward_, reverse_


def get_forward(wildcards):
    """Get forward reads for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward"].tolist()[0]


def get_reverse(wildcards):
    """Get reverse reads for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse"].tolist()[0]
