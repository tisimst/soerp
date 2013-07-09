from distutils.core import setup

with open('README.txt') as file:
    long_description = file.read()

setup(name='soerp',
    version='0.8.1',
    author='Abraham Lee',
    description='Second Order ERror Propagation',
    author_email='tisimst@gmail.com',
    url='http://pypi.python.org/pypi/soerp',
    license='BSD License',
    long_description=long_description,
    package_dir={'soerp':''},
    packages=['soerp'],
    requires=['ad'],
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
