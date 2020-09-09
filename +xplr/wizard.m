function wizard(file)

if nargin < 1
    file = brick.getfile('*.csv', 'Select data files to XPLOR');
end
data = io.read_table(file);
xplor(data)
