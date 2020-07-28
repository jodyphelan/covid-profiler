from setuptools import find_packages, setup

setup(
    name='covid-profiler',
    version='0.2.0',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'flask',
    ],
    scripts= [
        "scripts/covid-profiler.py",
        "scripts/covid_profiler_mask_fasta.py",
        "scripts/covid_profiler_mask_fasta_non_acgt.py"
    ],
)
