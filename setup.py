from setuptools import setup

setup(
    name='simobj',
    version='0.1',
    description='Code abstraction of objects (galaxies) in simulations.',
    url='',
    author='Kyle Oman',
    author_email='koman@astro.rug.nl',
    license='',
    packages=['simobj'],
    install_requires=['numpy', 'astropy', 'simfiles', 'utilities'],
    include_package_data=True,
    zip_safe=False
)
