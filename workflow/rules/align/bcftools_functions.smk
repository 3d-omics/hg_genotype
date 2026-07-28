BCFTOOLS_WINDOW_SIZE = 100_000
BCFTOOLS_WINDOW_PADDING = 1_000


def _chop_region_into_padded_windows(chrom, start, end):
    """Chop chrom:start-end into ~100kb windows, each padded on both sides.

    mpileup/call run over the padded interval so indel calling near a
    window edge has full read context, but only the core (unpadded)
    interval is kept in that window's output (`bcftools filter --targets
    ... --targets-overlap 0`; --targets streams rather than index-jumping,
    so it works on piped input). Adjacent windows' core intervals partition
    the region exactly, so every site is assigned to exactly one window and
    concatenation never sees duplicate or missing boundary sites.
    """
    windows = []
    win_start = start
    while win_start <= end:
        win_end = min(win_start + BCFTOOLS_WINDOW_SIZE - 1, end)
        pad_start = max(start, win_start - BCFTOOLS_WINDOW_PADDING)
        pad_end = min(end, win_end + BCFTOOLS_WINDOW_PADDING)
        windows.append(
            {
                "padded": f"{chrom}:{pad_start}-{pad_end}",
                "core": f"{chrom}:{win_start}-{win_end}",
            }
        )
        win_start = win_end + 1
    return windows


def _build_bcftools_windows():
    """Precompute, from REGIONS_BED4, one padded/core interval pair per
    ~100kb window of every region. Region lengths are already known at
    parse time (loaded straight from the BED4 file), so the whole window
    list can be enumerated here instead of behind a checkpoint.
    """
    windows = {}
    for region in REGIONS:
        row = REGIONS_BED4[REGIONS_BED4.name == region].iloc[0]
        chrom, start, end = row.chrom, int(row.chromStart), int(row.chromEnd)
        for i, window in enumerate(_chop_region_into_padded_windows(chrom, start, end)):
            windows[f"{region}__{i:05d}"] = window
    return windows


BCFTOOLS_WINDOWS = _build_bcftools_windows()
BCFTOOLS_WINDOW_NAMES = list(BCFTOOLS_WINDOWS)


def get_bcftools_window_padded(wildcards):
    return BCFTOOLS_WINDOWS[wildcards.window]["padded"]


def get_bcftools_window_core(wildcards):
    return BCFTOOLS_WINDOWS[wildcards.window]["core"]
