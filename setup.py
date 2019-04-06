from setuptools import setup

setup(
    name='ecp',
    url='https://github.com/jladan/package_demo',
    author='Ethan Goolish',
    author_email='ethangoolish@gmail.com',
    packages=['ecp'],
    install_requires=['numpy', 'scipy'],
    version='0.1',
    license='MIT',
    description='ecp is a Python package containing various procedures for finding multiple change-points',
    long_description=
    '''
    ecp is the Python equivalent of the ecp R package. It contains various procedures for finding multiple change-points.
    Two methods make use of dynamic programming and pruning, with no distributional assumptions other than the existence
    of certain absolute moments in one method. Hierarchical and exact search methods are included. All methods return the
    set of estimated change-points as well as other summary information.

    Please see package manual for more details.
    '''
)