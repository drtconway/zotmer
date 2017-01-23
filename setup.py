from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        pytest.main(self.test_args)

setup(name='zotmer',
        version='0.1',
        description='A python k-mer application workbench',
        url='http://github.com/drtconway/zotmer',
        author='Tom Conway',
        author_email='drtomc@gmail.com',
        license='Apache2',
        keywords='bioinformatics genomics pathogenomics',
        packages=find_packages(),
        install_requires=['docopt', 'pykmer>=0.1'],
        entry_points = {
            'console_scripts' : ['zot=zotmer.cli:main']
        },
        tests_require=['pytest'],
        cmdclass = {'test': PyTest},
        zip_safe=False)
