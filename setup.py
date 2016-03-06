from setuptools import setup

setup(name='pyfrp',
      version='1.0',
      description='PyFRAP: A Python based FRAP analysis tool box',
      url='https://github.com/alexblaessle/PyFRAP',
      author='Alexander Blaessle',
      author_email='alexander.blaessle@tuebingen.mpg.de',
      license='GNU GPL Version 3',
      packages=['pyfrp','pyfrp.modules','pyfrp.subclasses','pyfrp.gui'],
      package_dir={'pyfrp': 'pyfrp',
		   'pyfrp.modules': 'pyfrp/modules',
		   'pyfrp.gui' : 'pyfrp/gui'
		   },
      classifiers= [
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Biophysics/FRAP :: Analysis/Visualization',
        'Intended Audience :: Science/Research'
        ],
      platforms=['ALL'],
      keywords=["FRAP", "fluorescence",'recovery','after','photobleaching','reaction-diffusion','fitting'
              ],
      zip_safe=False)