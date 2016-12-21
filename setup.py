from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

setup(name='zotmer',
        version='0.1',
        description='A python k-mer application workbench',
        url='http://github.com/drtconway/pykmer',
        author='Tom Conway',
        author_email='drtomc@gmail.com',
        license='Apache2',
        keywords='bioinformatics genomics pathogenomics',
        packages=find_packages(),
        install_requires=['docopt', 'pykmer>=0.1'],
        entry_points = {
            'console_scripts' : ['zot=zotmer.cli:main']
        },
        zip_safe=False)
