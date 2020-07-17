function x = readmovie(video_file)

% file name
if nargin < 1
    video_file = brick.getfile;
end

video = VideoReader(video_file);
v = video.read();
rate = video.FrameRate;

% remove color
v = squeeze(v(:,:,1,:));


x = xplr.XData(v, ...
    {{'x' 'px' 1}, {'y' 'px' 1}, {'time' 's' 1/rate}});
    
    