def get_sex(sample):
    """From a sample name, get it's sex from the samples table"""
    sex = (
        samples[samples["sample"] == sample]["sex"].drop_duplicates().values.tolist()[0]
    )
    return sex


def get_sample_chromosomes(sample):
    """Get the chromosomes for a given sample"""
    sex = get_sex(sample)
    autosomes = features["reference"]["autosomes"]
    male_chromosomes = features["reference"]["male_chromosomes"]
    female_chromosomes = features["reference"]["female_chromosomes"]
    mitochondria = features["reference"]["mitochondria"]
    if sex == "male":
        return autosomes + male_chromosomes + mitochondria
    else:
        return autosomes + female_chromosomes + mitochondria
