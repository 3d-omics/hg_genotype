def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]
