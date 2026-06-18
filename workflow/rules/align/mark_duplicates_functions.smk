def _get_libraries_for_sample(sample_id):
    return samples[samples["sample_id"] == sample_id].library_id.values.tolist()


def get_crams_for_mark_duplicates(wildcards):
    """Get all CRAMs for a sample across libraries for mark_duplicates."""
    return [
        MAP / f"{wildcards.sample_id}.{library_id}.cram"
        for library_id in _get_libraries_for_sample(wildcards.sample_id)
    ]


def get_crais_for_mark_duplicates(wildcards):
    return [
        MAP / f"{wildcards.sample_id}.{library_id}.cram.crai"
        for library_id in _get_libraries_for_sample(wildcards.sample_id)
    ]
