
# xplor classes

## Main window

### `view`

**xplor** allows visualizing data of any dimensionality in windows that are linked to each other, using a syntax as simple as `xplor data`.

Each window corresponds to a `view` object, which is therefore the central object of the toolbox.
A `view` object contains data and operands (in other words, a `slicer`) and displays the resulting `slice` in multiple possible ways. Its main components are a panel (object `viewcontrol` below) to control data operations, and a main display (`viewdisplay`) for visualization.

Below we describe successively the classes used for data and header representation, operations, and different display 
components. Note that the **xplor** toolbox uses classes and functions from the **brick** toolbox, but we do not detail those here.


## Basic data and header


### `dimensionlabel`
A `dimensionlabel` is simply, as its name indicates, the label for a given data dimension. In addition to the label itself, it also indicates the class of values ('numeric', 'char' or 'logical'), and, in the case of numeric class and if appropriate, a unit of reference (e.g. 's') and a table of values for other possible units (e.g. 'ms': 1e-3, 'min': 60, 'hour': 3600).


### `header`
header information for a single data dimension; it entails a label (property `label`), the length of the data and values corresponding to each data point.

* In the general case, called _categorical_, values consist of a table (property `values`) with as may rows as data points. Also, each column is described by a `dimensionlabel` object stored in the `sublabels` property (whose length is therefore equal to the number of columns in the values table). Table of values and sublabels can also be left empty.
* In the particular case where data points are sampled along a regular grid (typically in time or space), the header is called a _measure_. Values are simply described by the start and step values. The `sublabel` property is a single numerical `dimensionlabel` object providing units information, and whose label should be identical to the header's own label.

`header` objects can't be modified, they can only be re-created  


### `xdata`
N-dimensional data and the corresponding N headers

Contrary to `header`, an `xdata` object is _active_ in the sense it can be modified, and sends notification upon modification, detailing what has changed (which dimension, which indices, etc.).


## Data operation and filtering


### `dataoperand (abstract)`
A `dataoperand` object describes an operation transforming some input data belonging to a given input space into an output data belonging to a possibly different output space.

Its `headerin` and `headerout` properties characterize these input and output spaces.

The description of its operation can be changed, in which case a particular event is generated.

The `dataoperand` cannot be instantiated (it is an abstract class), only subclasses defining specific operations can (see below).


### `point < dataoperand`
A `point` object is the simplest possible `dataoperand`, as it simply slices the data at a given position. The number of dimensions of its input and output are respectively 1 and 0 (i.e. output is a scalar).


### `filter < dataoperand`
A `filter` object is performing averaging in a given set of dimensions. Its `selections` property describes a set of regions of interets on which to perform this averaging. The number of regions of interest will determine the length of the output.

Filters applying to categorical data ressemble those in classical data visualization software (e.g. Excell, Tableau). Note that even though they apply on a single dimension, this dimension itself can be multi-label, therefore they in fact can filter several variables. [_not implemented yet_]


### `filterAndPoint < dataoperand`
The `filterAndPoint` is a combination of `point` and `filter` which processes data in the following way: if the set of selections is non-empty, it averages data over these ROIs, otherwise it uses the set of points to simply perform slicing.

### `zoomfilter < dataoperand`

_to be completed_

### `slicer`
A `slicer` object contains an `xdata` object and a set of `dataoperand` objects to be applied to specific dimensions of the data. Its output is stored in the `slice` property. It automatically handles notifications of changes of data or operands, as well as adding/removing operands, by modifying the `slice` accordingly.

Only operands that consist in averaging in orthogonal dimensions are accepted by the `slicer` object: the order in which such operations are performed does not matter,  therefore the `slicer` chooses the order that is the most efficient.

### `zoomslicer < slicer`

_to be completed, tell about the very tricky parts_

## Bank

### `bank`
The `bank` class has only one instantiation, as it repertoriates all allready-created `dimensionlabel`, `header`, and `dataoperand` objects, so as to ease linking of different data together.

It also remembers previous definitions from previous Matlab sessions, to avoid re-defining such objects. 


## Controls

Graphic objects can be coarsely divided into two groups, depending on whether they are more dedicated to control what is being displayed and how, or to actually display data. Of course "control" objects do display some information already, and "display" objects also exert control on other linked displays.

### `viewcontrol`

The `viewcontrol` object is the part of the `view` window that controls how the data is displayed. It consists in a panel where one can choose existing or new filters to be applied to each dimension (or to a set of dimensions, e.g. when creating 2D ROIs). New created filters will be either  _shared_ or _private_: _shared_ filters are visible by other `view` objects and their controls (typically a `list` object) will be opened in a new window. _Private_ filters apply only to the data inside the considered `view` window and are opened inside a sub-panel of the `viewcontrol` panel.


### `list`

A `list` is as its name indicates a list, acting on an unidimensional `filterAndPoint` filter.


Multiple selections can be made and are marked with numbers under brackets inside the list display. A drop-down menu (activated by right-clicking the list display) offers multiple options for making single-point or multiple-points selections. An important distinction exists between _solid_ selections made through this menu and _soft_ selections made by simply selecting entries in the list: contrary to the solid ones, these _soft_ selections are removed as soon as a new selection is made.

_Under development: a write-only property `singleSelection` will specify whether only one selection can be made. In that case the `list` behavior will be simpler since the unique selection will be simply defined by the selected elements in the list._

### listcombo

A `listcombo` is a window or panel grouping several `list` controls.

## Display

### `viewdisplay`

The `viewdisplay` object displays data into a graph. Currently available modes are 'time courses' and 'image', which can represent arbitrary N-dimensional data by tiling (for example a 6D array can be visualized as a 2D array of 2D arrays of images). The `org` property describes how this visualization is arranged.

Separate classes (see below) take charge of specific objects inside the display, and communicate with the parent `viewdisplay` through method calls.


### `displaylabels`  

The `displaylabels` object automatically updates display of data header labels upon changes of the parent `viewdisplay` handling of different dimensions (`org` and `activedim` properties), and conversely offers user mouse control on these properties.

Each label has its own context menu that gives user control over properties that are specific to each dimension (e.g. binning value, etc.). Therefore the `displaylabels` class also controls such properties. To do so, it relies on another class, `displaylabelmenu`: a `displaylabels` object contains one `displaylabelmenu` object per label, i.e. per data dimension.


### `displaygraph`

The `displaygraph` object controls the organization of the axes where data will be displayed. In particular it offers conversion facilities between N-dimensional coordinates of data points and 2-dimensional coordinates of positions in the axes where these data points are displayed. It also takes care of maintaining appropriate axis ratio when needed (e.g. to make image pixels appear square). 