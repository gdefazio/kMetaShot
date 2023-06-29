from setuptools import setup, find_packages

setup(
    name='kmetashot',
    version='1.0',
    scripts=['/home/giuseppedefazio/Documenti/tesi_dottorale/kMetaShot/kMetaShot_classifier_NV.py',
             '/home/giuseppedefazio/Documenti/tesi_dottorale/kMetaShot/test/kMetaShot_test.py'],
    packages=['kMetaShot_package'],
    install_requires=['numpy==1.18.1',
                      'numba==0.51.2',
                      'wget',
                      'pandas==1.0.4',
                      'mmh3',
                      'bitarray===1.2.1',
                      'h5py==2.9.0',
                      #'hdf5==1.10.4',
                      'argcomplete==1.11.1'],
    url='https://github.com/gdefazio/kMetaShot',
    license='GNU GENERAL PUBLIC LICENSE',
    author='Defazio G. et al.',
    author_email='bruno.fosso@uniba.it',
    description='''
    kMetaShot is a taxonomy classifier tool for 2nd and 3rd GS data for 
    shotgun metagenomics relying on an aligmnent free and k-mer/minimizer 
    based approach.
    '''
)
