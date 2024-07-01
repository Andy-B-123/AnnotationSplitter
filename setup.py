# setup.py

from setuptools import setup, find_packages

setup(
    name="AnnotationSplitter",
    version="0.1.0",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        "biopython>1",
        "cogent3>=2024.4",
        "pandas>1",
        "matplotlib>3",
        "seaborn>0.10"
    ],
    entry_points={
        'console_scripts': [
            'AnnotationSplitter=main:main',
        ],
    },
    author="Andreas Bachler",
    author_email="Andy.Bachler@example.com",
    description="A simple bioinformatics script.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/Andy-B-123/AnnotationSplitter",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.11',
)
