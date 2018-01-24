# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.txt') as f:
    readme = f.read()

with open('LICENSE.txt') as f:
    license = f.read()

setup(
    name='des_stacks',
    version='0.1.0',
    description='Framework for deep coadds of DES SN fields',
    long_description=readme,
    author='Philip Wiseman',
    author_email='pacelweb@gmail.com',
    url='https://github.com/wisemanp/des_stacks',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)
