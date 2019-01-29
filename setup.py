import setuptools
from setuptools_rust import Binding, RustExtension

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='alignmentrs',
    author='Kent Kawashima',
    version='0.1.0',
    author_email='kentkawashima@gmail.com',
    description='Quickly read and manipulate multiple sequence alignments in Python',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/kentwait/alignmentrs',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords=['block', 'alignment', 'bioinformatics'],
    rust_extensions=[
        RustExtension('libalignmentrs.sequence',
                      'Cargo.toml', binding=Binding.PyO3),
        RustExtension('libalignmentrs.marker',
                      'Cargo.toml', binding=Binding.PyO3),
    ],
    packages=['alignmentrs'],
    package_data={
        'alignmentrs': ['lib/libalignmentrs/alignment.cpython-37m-darwin.so']},
    requires=[],
    zip_safe=False,  # Rust extensions are not zip safe, like C-extensions.
)
