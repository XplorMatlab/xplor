classdef zoomcentral < xplr.graphnode
    % zoom central is made to synchronized zoom, it needs to be an external
    % object to use as reference to sync different zooms with differents
    % scale of the same unit
    %   zoomcentral herit from zoomfilter,
    
    properties (SetAccess='private')
        zoom
        label
        unit
    end
    
    events
       ChangedZoom
    end
    
    methods
        function ZC = zoomcentral(label,unit)
            % zoomcentral(label,unit)
            % To create a zoomcentral
            ZC.label = label;
            ZC.unit = unit;
        end
        
        function connectZoomFilter(ZC,zoomFilterToConnect)
            % connect a new zoomfilter to the zoomcentral using graphnode
            % listeners connector
            ZC.addListenerExclusivePair(zoomFilterToConnect, ...
                'ChangedZoom',@(u,e)ZC.centralToLocal(zoomFilterToConnect), ...
                'ChangedZoom',@(u,e)ZC.localToCentral(zoomFilterToConnect));
        end
        
        function localToCentral(ZC,zoomFilter)
            %get scale
            zoomFilterScale = zoomFilter.headerin.scale;
            zoomFilterStart = zoomFilter.headerin.start;
            
            % TODO: bug fix, this is a special case when the zoom is ":".
            % The correct behavior has to be defined when zoom is reset
            % 
            if(zoomFilter.zoom == ':')
                zoomToSet = ':';
            else
                zoomToSet = zoomFilterStart + (zoomFilter.zoom-1)*zoomFilterScale;
            end
            %sourceScale = ZC.label.scale;

            ZC.setZoom(zoomToSet);
        end
        
          % synchronize zoom with scale adaptation
        function centralToLocal(ZC,zoomFilter)
            zoomFilterScale = zoomFilter.headerin.scale;
            zoomFilterStart = zoomFilter.headerin.start;
            
            % TODO: bug fix, this is a special case when the zoom is ":".
            % The correct behavior has to be defined when zoom is reset
            % 
            if(ZC.zoom == ':')
                zoomToSet = ':';
            else
                zoomToSet = 1 + (ZC.zoom - zoomFilterStart)/zoomFilterScale;
            end

            zoomFilter.setZoom(zoomToSet);
        end
        
        function setZoom(ZC,zoom)
            ZC.zoom = zoom;
            notify(ZC,'ChangedZoom');
        end
    end
end