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


GENOTYPE = Path("results/genotype/")
GATK = GENOTYPE / "gatk"

ANNOTATE = Path("results/annotate/")
SNPEFF = ANNOTATE / "snpeff"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
REPORT_CHR = Path("reports/by_chromosome/")

SWAPS = Path("results/swaps/")
