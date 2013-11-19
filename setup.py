import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='soerp',
    version='0.9.5',
    author='Abraham Lee',
    description='Second Order Error Propagation',
    author_email='tisimst@gmail.com',
    url='https://github.com/tisimst/soerp',
    license='BSD License',
    long_description=read('README.rst'),
    packages=[
        'soerp', 
        'soerp.umath'
        ],
    install_requires=[
        'ad', 
        'numpy',
        'scipy',
        'matplotlib',
        ],
    keywords=[
        'uncertainty analysis', 
        'uncertainties', 
        'error propagation', 
        'second order', 
        'derivative', 
        'statistics', 
        'method of moments', 
        'distribution'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.0',
        'Programming Language :: Python :: 3.1',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Topic :: Education',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development',
        'Topic :: Software Development :: Libraries',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: Utilities'
        ]
    )
