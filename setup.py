import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scraper",
    version="0.0.1",
    author="Katarina Gresova",
    author_email="gresova11@gmail.com",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/katarinagresova/ensembl_scraper",
    packages=setuptools.find_packages(include=['scraper']),
    python_requires='>=3.6',
)