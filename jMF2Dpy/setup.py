import setuptools

__lib_name__ = "jMF2D"
__lib_version__ = "1.1.0"
__description__ = "jMF2D: Enhancing and accelerating cell type deconvolution " \
                  "of spatial transcriptomics with dual network model"
__url__ = "https://github.com/xkmaxidian/jMF2D"
__author__ = "Yuhong Zha"
__author_email__ = "18583330027@163.com"
__license__ = "MIT"
__keywords__ = ["Spatial transcriptomics", "Deconvolution", "Joint learning", "Topological structure", "Integrative analysis"]
__requires__ = [
    "pandas==2.0.3",
    "numpy==1.24.4",
    "scanpy==1.9.8",
    "scipy==1.10.1",
    "scikit-learn==1.0.2",
    "scikit-misc==0.2.0",
    "matplotlib==3.7.5",
    "seaborn==0.13.2",
    "anndata==0.9.2",
    "opencv-python==4.10.0.82",
    "leidenalg==0.10.2",
    "POT==0.9.3",
]
setuptools.setup(
    name=__lib_name__,
    version=__lib_version__,
    author=__author__,
    author_email=__author_email__,
    description=__description__,
    long_description_content_type="text/markdown",
    url=__url__,
    install_requires=__requires__,
    packages=setuptools.find_packages(),
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)