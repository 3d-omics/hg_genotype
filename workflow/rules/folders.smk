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
SNPEFF_DB = Path("resources/snpeff/")
SNPEFF = ANNOTATE / "snpeff"


REPORT = Path("reports/")
STEP = REPORT / "step"
LIBRARY = REPORT / "library"

SWAPS = Path("results/swaps/")
