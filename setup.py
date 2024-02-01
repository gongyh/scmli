try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup


setup(
    name="scmli",
    version="0.0.1",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    packages=find_packages(),
    python_requires=">=3",
    install_requires=[
        "biopython",
        "pandas",
        "lxml",
    ],
    entry_points={"console_scripts": ["scmli = scmli:main"]},
)
