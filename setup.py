import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FLAME: Full Length Adjecency Matrix Enumeration", # Replace with your own username
    version="0.3.1",
    author="Alan Bäckerholm",
    author_email="Alan.baek1@gmail.com",
    description="FLAME: Full Length Adjecency Matrix Enumeration - is a longread splice variant annotation tool that allows for the filtering and quantification of splice variants according to the gene annotation file. It also allows for the detection and confirmation of novel exons and/or splice variants through the use of adjecent nucleotide in the reference and through the use of shortread sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
