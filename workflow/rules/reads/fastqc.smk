rule reads__fastqc__all:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample_id}.{library_id}_{end}_fastqc.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],
