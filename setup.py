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
    package_data={'talos_af': ['example_config.toml']},
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
    include_package_data=True,
    install_requires=[],
    extras_require={
        'test': [],
    },
    entry_points={
        'console_scripts': [
        ],
    },
)