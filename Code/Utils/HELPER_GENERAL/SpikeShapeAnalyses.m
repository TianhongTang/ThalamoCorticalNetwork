
if PDS.plexonconv(1)==0
    
    try
        load(['mapfile' mat2str(xlsraw{xyz,13}) '.mat']); %get spike shape
        wave = CSPK_001_BitResolution*  mean(CSEG_001_Template1_SEG');
        wave1=wave;
        ymiN = find( wave1 == min(wave1) ); ymiN=ymiN(1);
        wave1(1:ymiN)=NaN;
        ypeaK = find( wave1 == max(wave1) ); ypeaK=ypeaK(1);
        wduRR= (ypeaK(1)-ymiN(1));
        wave=wave-min(wave)
        wave=wave./max(wave)
        
        
        savestruct(xyz).wave=wave(1:92);
        savestruct(xyz).waveduration= wduRR./44;
        savestruct(xyz).ypeaK= ypeaK;
        savestruct(xyz).ymiN= ymiN;
        
    catch
        
        savestruct(xyz).wave=NaN;
        savestruct(xyz).waveduration= NaN;
        savestruct(xyz).ypeaK= NaN;
        savestruct(xyz).ymiN= NaN;
        
    end
    
else
    
    try
        
        wave=PDS.waveshape;
        wave1=PDS.waveshape;
        ymiN = find( wave1 == min(wave1) ); ymiN=ymiN(1);
        wave1(1:ymiN)=NaN;
        ypeaK = find( wave1 == max(wave1) ); ypeaK=ypeaK(1);
        wduRR= (ypeaK(1)-ymiN(1));
        wave=wave-min(wave)
        wave=wave./max(wave)
        wave(length(wave)+1:92)=NaN;
        
        savestruct(xyz).wave=wave(1:92);
        savestruct(xyz).waveduration= wduRR./PDS.waveshapeSMP;
        savestruct(xyz).ypeaK= ypeaK;
        savestruct(xyz).ymiN= ymiN;
        
    catch
        
        savestruct(xyz).wave=NaN;
        savestruct(xyz).waveduration= NaN;
        savestruct(xyz).ypeaK=NaN;
        savestruct(xyz).ymiN= NaN;
        
    end
    
    
end


clear spikeshape wave wave1 wduRR CSPK_001_BitResolution CSEG_001_Template1_SEG