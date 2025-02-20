from setuptools import setup, find_packages

setup(
    name='gosh',
    version='0.1.0',
    description='gOSh - gOS sHell: A CLI tool for Nextflow pipeline management',
    author='MSK Lab',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'openai',
    ],
    entry_points={
        'console_scripts': [
            'gosh=gosh.main:cli',
        ],
    },
    python_requires='>=3.6',
)
