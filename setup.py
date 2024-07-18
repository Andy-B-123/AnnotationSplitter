from setuptools import setup, find_packages

setup(
    name="AnnotationSplitter",
    version="0.0.4",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        "icecream>=2.1",
        "pandas",
        "tqdm",
        "requests",
        "cogent3>=2024.2.5a1",
        "pyfaidx>=0.8",
        "biopython>=1.83"
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
