
from distutils.core import setup

setup(name='pyLMM',
      version='0.99b',
      author="Nick Furlotte",
      author_email='nick.furlotte@gmail.com',
      url="https://github.com/nickFurlotte/pylmm",
      description='pyLMM is a lightweight linear mixed model solver for use in GWAS.',
      packages=['pylmm'],
      scripts=['scripts/pylmmGWAS.py', 'scripts/pylmmKinship.py',  'scripts/estimate_variance_components.py'],
      )
