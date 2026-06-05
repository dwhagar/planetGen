import nltk
import os
from setuptools import setup, find_packages
from setuptools.command.install import install


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        # Run the standard install
        install.run(self)
        # Download the NLTK 'words' corpus
        print("Downloading NLTK 'words' corpus...")
        nltk.download('words', quiet=True)
        print("NLTK 'words' corpus downloaded successfully.")


setup(
    name='planetGen',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'nltk',
        'transformers',
        'torch'
    ],
    cmdclass={
        'install': PostInstallCommand,
    },
    author='David Hagar',
    author_email='david.hagar@gmail.com',
    description='A procedural planet and star system generator.',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/planetGen',  # Replace with your project's URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
