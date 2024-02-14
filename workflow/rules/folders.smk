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
POSTERIORS = GENOTYPE / "variant_posteriors"
COMBINE_GVCFS = GENOTYPE / "combine_gvcfs"
HAPLOTYPE_CALLER = GENOTYPE / "haplotype_caller"
GENOTYPE_GVCFS = GENOTYPE / "genotype_gvcfs"
VARIANT_FILTRATION = GENOTYPE / "variant_filtration"

ANNOTATE = Path("results/annotate/")
SNPEFF_DB = Path("resources/snpeff/")
SNPEFF = ANNOTATE / "snpeff"


REPORT = Path("reports/")
STEP = REPORT / "step"
LIBRARY = REPORT / "library"

SWAPS = Path("results/swaps/")
