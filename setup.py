from setuptools import setup, find_packages

print(find_packages(where='src'))  # Add this line to check the output

setup(
    name="AnnotationSplitter",
    version="0.1.6",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        "scipy>=1.13.1",
        "pandas>=2.2.2",
        "numpy>=1.26.4",
        "matplotlib>=3.9.0",
        "cogent3>=2024.5.7a1",
        "tqdm>=4.66.4",
        "icecream>=2.1.3",
        "biopython>=1.8",
        "pyfaidx>=0.8.1.1",
        "requests",
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
