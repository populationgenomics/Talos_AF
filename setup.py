"""
setup.py for the talos package
"""

from setuptools import find_packages, setup

with open('README.md', encoding='utf-8') as handle:
    readme = handle.read()


setup(
    name='talos_af',
    description='Centre for Population Genomics: Actionable Findings Detection',
    long_description=readme,
    version='0.0.0',
    author='Matthew Welland, CPG',
    author_email='matthew.welland@populationgenomics.org.au, cas.simons@populationgenomics.org.au',
    # include_package_data=True,
    # package_data={'talos_af': ['example_config.toml']},
    url='https://github.com/populationgenomics/talos_af',
    license='MIT',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        'bump2version>=1.0.1',
        'cloudpathlib[all]>=0.18.1',
        'cyvcf2==0.31.1',
        'hail==0.2.133',
        'httpx==0.27.0',
        'peds==1.3.2',
        'pre-commit>=3.7',
        'pydantic==2.5.2',
        'ruff>=0.4.7',
        'tenacity>=9.0.0',
        'toml==0.10.2',
    ],
    extras_require={
        'test': [
            'pytest>=8.3',
        ],
    },
    entry_points={
        'console_scripts': [],
    },
)
