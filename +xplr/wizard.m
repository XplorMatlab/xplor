
file = brick.getfile('*.csv', 'Select data files to XPLOR');
data = io.readtable(file);
xplor(data)
