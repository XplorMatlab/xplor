loading the data
++++++++++++++++

The first thing to do is to load your data. Your data must be stored into a N
dimension array.

Then you can launch the toolbox with this simple line of command : xplor mydata

**Screenshot of the command line to start xplor**

	.. image:: image/callingxplor.png
	  :width: 400
	  :alt: Set Path

The following window appears. It allows you to enter the metadata. For each
dimension, there is a label, possibly a unit.
The scale/value column can be filled with :

- nothing (the ticks of the axis will simply be numbered from 0)

- a list of named elements (each name will appear with a tick on the axis)

- a start and a scale number (the ticks will be calculated from those
   numbers), it corresponds to a continous dimension.

setting the headers
+++++++++++++++++++

**Screenshot of the set headers window**

	.. image:: image/setheaders.png
	  :width: 600
	  :alt: Set Path

Some filters may be automatically applied.
This will reduce the size of the data, which will then be displayed properly.

**Screenshot of the initial display**

	.. image:: image/initialdisplay.png
	  :width: 600
	  :alt: Set Path

create a new filter
+++++++++++++++++++

To create a filter, right click on the name of the dimension you wish to
filter in the left part of the window. A menu will appear. Click on Filter with share 1D filter.

Each time a filter is created, a list of all the elements of the dimension
(what appears on the axis) appears in the Shared Filters window.
This window allows you to choose what elements to display on the graph
(see next paragraphs).
All the filters also appear on the main window in purple. It informs you of the
order in which the filters are applied, and allows you to change it as well as
locally deactivate a filter by clicking on it.


**Screenshot of the creation of filters**

	.. image:: image/filtercreation.png
	  :width: 800
	  :alt: Set Path

In the shared filter window, you can select for each filtered dimension the
values that will be displayed.

For instance, in the previous screenshot we had the unemployment rate for
females aged between 15 and 19 years old.

By clicking on F and M in the sex column and Y15-39 and Y40-64 in the age
column, one can easily compare the trends as seen in the next display.


Change the value of a filter
++++++++++++++++++++++++++++

**Screenshot of a comparative display**

	.. image:: image/comparison.png
	  :width: 600
	  :alt: Set Path

If you which to display the same imformation differently, you can drag the
labels. On the next screenshot, "sex" is not on the same axis as before.


Average some values
+++++++++++++++++++

**Screenshot of an average**

	.. image:: image/average.png
	  :width: 600
	  :alt: Set Path

To make an average, as on the previous screenshot, select the values you want
to average. Right click on it and a menu will appear. Click on "new group
selection". The values are now grouped, you can then compare them to another
single selection or another group selection.

