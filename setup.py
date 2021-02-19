#!/usr/bin/env python
# coding: utf-8

from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    readme = fh.read()

setup(
    name='h1d',
    version='0.0.25',
    author='wangjiankng',
    author_email='wangjk321@gmail.com',
    url='https://github.com/wangjk321/HiC1Dmetrics',
    description='HiC1Dmetrics pip version',
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=['h1d'],
    include_package_data=True,
    install_requires=["pandas","numpy","scikit-learn",
		"fithic==2.0.7","statsmodels"],
    classifiers=[
        "Programming Language :: Python :: 3",
	"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    entry_points={
        'console_scripts': [
            'h1d=h1d.__main__:CLI',
        ]}
)
