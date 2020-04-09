.. toctree::
   :maxdepth: 2
   :caption: Contents:


Structure
*********

The three main modules are :
   - **view**, the user interface,
   - **xdata**, the data container,
   - **operation**, to apply filters.


The following UML diagram shows the the main structure of the toolbox. The use
of each of the module is presented below. They are further detailed in the
next pages of the documentation.
At the last section of the page, the main interactions between the classes are
described in sequence diagrams.

UML diagram of all Xplor's modules
----------------------------------

**UML diagram of the global structure of the toolbox**

.. image:: \images\UML\GlobalDiagram.png
   :align: center
   :alt: UML diagram of all Xplor's modules


Presentation of the modules
---------------------------

view
++++

**view**'s main class is Window. Window creates the canvas (i.e. the window)
on witch to display the data. This canvas is composed of a control zone, where
the user can select the filters (of the operation module) to apply, and a
display zone, in which the user can zoom and change the labels' position.
Window also possesses a Slicer instance and a ZoomSlicer instance (both of the
operation module) to apply the correct filters on the Xdata instance.

**Screenshot of the Window with the control zone and the display zone**

.. image:: \images\screenshots\basicexample\windowExp.png
   :align: center
   :alt: Screenshot of Window

**Interactions between the Window, the Slicer, the Filters and the Xdata**

.. image:: \images\drawings\slicerandfilter.png
   :align: center
   :alt: Main mecanism of Xplor

xdata
+++++

The module **xdata** aims at the creation of a container for the N dimensional
data and all the relative information on the dimensions and the data itself.

The metadata of a dimension is stored in a header.

**Example of a MeasureHeader instance**

.. image:: \images\examples_module_xdata\measureheader.png
   :align: center
   :alt: illustration for a measure header

The data, the headers and other metadata are contained in the xdata object.

**Example of a Xdata instance**

.. image:: \images\examples_module_xdata\xdataexample.png
   :align: center
   :alt: illustration for a xdata element

operation
+++++++++

The module **operation** contains two main classes :
   - Filters
   - Slicers

The slicer and zoomslicer elements apply a succession of Filters on
the Xdata instance. The result of those operations is what is going to be
displayed on the canvas.

The interactions between the various classes are illustrated by the sequence
diagrams of the next section.

list_display
++++++++++++

**list_display**'s main class is ListDisplay. It is an interface for the user
to select what data is to be displayed. For each filter created from the
control zone of a Window instance, a new List instance appears.

**Screenshot of the ListDisplay and the Window**

.. image:: \images\screenshots\basicexample\Filtre.png
   :align: center
   :alt: Screenshot of the ListDisplay element

bank
++++

**bank** is used to store previously used Headers and units, to save the user
some time.

header_edit
+++++++++++

**header_edit** allows the initialisation of a Xdata element : the user
specifies the headers' labels and units. It suggests some previously used
labels and units thanks to Bank.

**Screenshot of the set headers window**

.. image:: \images\screenshots\unemploymentexample\setHeaders.PNG
   :align: center
   :alt: Screenshot of set headers window


GUI_tools
+++++++++

**GUI_tools** is a set of tools for the interface.


Sequence Diagrams, events handling
----------------------------------

Initialisation of the Window
++++++++++++++++++++++++++++

**Sequence Diagram of the creation of a Window**

.. image:: \images\sequence\InitialisationSequence.png
   :align: center
   :alt: Sequence diagram of the creation of a Window

**Screenshot of the Window with the control zone and the display zone**

.. image:: \images\screenshots\basicexample\windowExp.png
   :align: center
   :alt: Screenshot of Window

In the control part of Window, a list of the dimensions appears. The user will
then be able to create a filter for each dimension by selecting it in the list.
What happens then is described in the next section.


Initialisation of a filter
++++++++++++++++++++++++++

**Sequence Diagram of the creation of a filter**

.. image:: \images\sequence\FilterSequence.png
   :align: center
   :alt: Sequence diagram of the creation of a filter

When the user creates a filter, it appears in the control zone as well as in
the ListDisplay element.

The ListDisplay element could interact with several windows. It's main function
is to allow the user to chose what values to give to a filter, to select what
data will be displayed.

The filter that appears in the control zone allows the user to change the
order of the filters to apply but also to activate/deactivate a filter for
this specific window.

Each time a filter is created, it's default value is to only keep the first
element of it's dimension.

**Screenshot of the newly created ListDisplay element and the Window**

.. image:: \images\screenshots\basicexample\ListDisplay.PNG
   :align: center
   :alt: Screenshot of ListDisplay and the Window

Changing a filter with the ListDisplay element
++++++++++++++++++++++++++++++++++++++++++++++

**Sequence Diagram of the modification of a filter**

.. image:: \images\sequence\ListDisplaySequence.png
   :align: center
   :alt: Screenshot of the ListDisplay

The user can change the selection in ListDisplay to change the filter. The new
selected slice will appear on the display zone.

**Screenshot of the changed ListDisplay element and the Window**

.. image:: \images\screenshots\basicexample\Filtre.PNG
   :align: center
   :alt: Screenshot of ListDisplay and the Window

Zooming in or out
+++++++++++++++++

**Sequence Diagram for a zoom action**

.. image:: \images\sequence\ZoomSequence.png
   :align: center
   :alt: Sequence diagram for a zoom action

**Screenshot of ListDisplay and the Window after a zoom action**

.. image:: \images\screenshots\basicexample\Zoom.PNG
   :align: center
   :alt: Screenshot after a zoom action

The side and upper bar allow the user to zoom in and out in each of the
dimensions separately. Only a part of the data will appear.

Another possibility is to bind the data in a dimension. All the elements
within the binding step will be averaged to form a new element. All the new
elements will be displayed.