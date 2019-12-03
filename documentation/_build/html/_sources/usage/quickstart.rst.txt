Quick Start
****************************

After :ref:`Install via Github`, you can load data into xplor:

* Use default datasets

	To try xplor, you can use the default datasets available in "xplor-Matlab-/demo data":
		- intrinsic.mat
		- unemployment.mat

	.. code-block:: javascript

		load('intrinsic.mat')
		xplor x

	or

	.. code-block:: javascript

		load('unemployment.mat')
		xplor dataall

* Use your own datasets

	To use your own datasets you first have to convert it to a N-dimensionnal array.
