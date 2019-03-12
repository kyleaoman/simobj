from setuptools import setup

setup(
    name='simobj',
    version='1.0',
    description='Code abstraction of objects (galaxies) in simulations.',
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='GNU GPL v3',
    packages=['simobj'],
    install_requires=['numpy', 'astropy', 'h5py', 'simfiles'],
    include_package_data=True,
    zip_safe=False
)
