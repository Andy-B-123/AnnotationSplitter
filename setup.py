from setuptools import setup, find_packages

setup(
    name="AnnotationSplitter",
    version="0.1.0",
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'AnnotationSplitter=src.main:main',
        ],
    },
    author="Andreas Bachler",
    author_email="Andy.Bachler@example.com",
    description="A simple Python project example.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/Andy-B-123/AnnotationSplitter",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
