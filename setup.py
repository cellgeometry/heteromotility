from setuptools import setup

setup(name='heteromotility',
      version='0.1.9',
      description='Quantitative analysis of biological motion',
      url='http://github.com/jacobkimmel/heteromotility',
      author='Jacob Kimmel',
      author_email='jacobkimmel@gmail.com',
      license='MIT',
      packages=['heteromotility'],
      install_requires = [
      'numpy',
      'scipy',
      'argparse',
      'statsmodels',
      'pandas',
      'patsy'
      ],
      entry_points = {
      'console_scripts': ['heteromotility = heteromotility.command_line:cli']
      },
      include_package_data = True,
      zip_safe=False)
