function coor = posPlot(yPos, xPos)
    width = 45;
    height = 25;
    
    coor.x1coor = width*(xPos-1) + 5;
    coor.x2coor = coor.x1coor + width-10;
    
    coor.y1coor = height*(yPos-1) + 5;
    coor.y2coor = coor.y1coor + height - 10;
end