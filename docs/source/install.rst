
Installation guide
==================

Install possibilities
---------------------
Below are several installation possibilities.

Install from PyPI (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To install the stable version from `PyPI <https://pypi.org/project/oscilate/>`_, use::

    pip install oscilate

Then, simply import the package in a Python environment using::

    import oscilate


Install from a GitHub release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To install from a GitHub release tagged as version `vX.Y.Z`, run::

    pip install https://github.com/vinceECN/OSCILATE/archive/refs/tags/vX.Y.Z.tar.gz



Install from the repository (latest version)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To install the latest version directly from the GitHub repository, run::

    git clone https://github.com/vinceECN/OSCILATE.git
    cd OSCILATE
    pip install .


Dependencies
------------

- **Python 3.8 or higher** is required.

- For development or building documentation, install additional dependencies::
  
    pip install -r requirements-dev.txt
    pip install -r docs/requirements.txt
  

Optional: use a virtual environment (recommended)
-------------------------------------------------
To avoid conflicts with other packages, create and activate a virtual environment::

    python -m venv venv_mms        ## Create the venv_mms virtual environment 
    source venv_mms/bin/activate   ## Linux/macOS
    .\venv_mms\Scripts\activate    ## Windows



Test the install
----------------

To test the install, follow these steps:

1. Open a Python environment. Ideally one powered by Jupyter (see :ref:`Outputs <outputs_section>` section) to display results as :math:`\LaTeX`.

2. In the documentation, go to :doc:`*Application Examples/Example 1* <examples/example1>`,

3. Copy the example code,

4. Run the example code in your Python environment.

5. You should see information about the ongoing computations.

6. After the code is ran (a few seconds should be enough), figures of the forced response and its stability information are displayed.

7. To access the analytical solutions computed, type, for instance, ``ss.sol.fa``. They will be displayed as :math:`\LaTeX` if the Python environment supports its. 
