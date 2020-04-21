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
        "scripts/update_tree.py",
        "scripts/list_mutations_by_node.py",
        "scripts/correct_covid_csq.py",
        "scripts/process-gisaid.py",
        "scripts/process_alignment.py"
    ],
)
