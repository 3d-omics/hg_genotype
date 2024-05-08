READS = Path("results/reads/")
REFERENCE = Path("results/reference/")

PRE = Path("results/preprocess/")
FASTP = PRE / "fastp"

ALIGN = Path("results/align/")
INDEX = ALIGN / "index"
MAP = ALIGN / "map"
SPLIT = ALIGN / "split"
MARK_DUPLICATES = ALIGN / "mark_duplicates"
RECALIBRATE = ALIGN / "recalibrate"
MERGE = ALIGN / "merge"


VARIANTS = Path("results/variants/")
CALL = VARIANTS / "call"
GENOTYPE = VARIANTS / "genotype"
FILTER = VARIANTS / "filter"
POSTERIORS = VARIANTS / "posteriors"

ANNOTATE = Path("results/annotate/")
SNPEFF_DB = Path("resources/snpeff/")
SNPEFF = ANNOTATE / "snpeff"
VEP = ANNOTATE / "vep"

REPORT = Path("reports/")
STEP = REPORT / "step"
SAMPLE = REPORT / "sample"

SWAPS = Path("results/swaps/")
SOMALIER = SWAPS / "somalier"
