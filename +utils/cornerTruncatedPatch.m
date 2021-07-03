function p = cornerTruncatedPatch(patchLength, patchWidth, cutLength, wavePolarization)
    if cutLength > max(patchLength, patchWidth)
        error("Cut too long");
    end
    
    if (wavePolarization ~= 'RHCP') & (wavePolarization ~= 'LHCP')
        error("wavePolarization must either be RHCP or LHCP");
    end
    
    p = antenna.Rectangle('Width', patchWidth, 'Length', patchLength);
    
    cornerVertices = [ patchLength/2,           -patchWidth/2 ;
                       patchLength/2-cutLength, -patchWidth/2 ; 
                       patchLength/2,           -patchWidth/2+cutLength] ;
    
    % Adjust polarization for LHCP
    if wavePolarization == 'LHCP'
        cornerVertices(:,1) = -cornerVertices(:,1);
    end
    
    p = p - antenna.Polygon('Vertices', cornerVertices);
    p = p - antenna.Polygon('Vertices', -cornerVertices);
end