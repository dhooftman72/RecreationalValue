%Phython codes for Recreational Value

% # 5km
% import arcpy
% arcpy.SetLogHistory(False)
% PointsShape = "5kmPoints" 
% SpeedRaster ="Speedraster_5Km.tif" 
% expression=" DH_ID > -9999"
% OutBase = "D:/Data/Dropbox/Lactuca/Temp/DistancePoint_"
% f = ["SelectedPoint","CostPoint","Divide_Cost"]
% for x in f:
%     arcpy.Delete_management(x)
% with arcpy.da.SearchCursor(PointsShape, ["DH_ID"], where_clause=expression) as cursor: 
%     for row in cursor:   
%           # Name and select Point
%           Point = '"DH_ID" = ' + str(row[0])  
%           with arcpy.EnvManager(overwriteOutput=True):
%               arcpy.analysis.Select(PointsShape, "SelectedPoint", Point)
%               out_distance_raster = arcpy.sa.CostDistance("SelectedPoint", "Speedraster_5Km.tif", None, None, None, None, None, None, '') 
%               out_distance_raster.save("CostPoint")
%               out_raster = arcpy.sa.Divide("CostPoint", 1000)
%               out_raster.save("Divide_Cost")
%               outfile = OutBase + str(row[0])+ ".asc"
%               arcpy.conversion.RasterToASCII("Divide_Cost", outfile)
%           f = ["SelectedPoint","CostPoint","Divide_Cost"]
%           for x in f:
%               arcpy.Delete_management(x)
%%
%%
% #2.5km with double expression
% import arcpy
% arcpy.SetLogHistory(False)
% PointsShape = "2_5kmPoints" 
% SpeedRaster ="Speedraster_2_5Km.tif" 
% expressionLow=" DH_ID > -9999 and  DH_ID < 50001"
% OutBase = "D:/Data/Dropbox/Lactuca/Temp/DistancePoint2_5_"
% f = ["SelectedPoint","CostPoint","Divide_Cost"]
% for x in f:
%     arcpy.Delete_management(x)
% with arcpy.da.SearchCursor(PointsShape, ["DH_ID"], where_clause=expressionLow) as cursor: 
%     for row in cursor:   
%           # Name and select Point
%           Point = '"DH_ID" = ' + str(row[0])  
%           with arcpy.EnvManager(overwriteOutput=True):
%               arcpy.analysis.Select(PointsShape, "SelectedPoint", Point)
%               out_distance_raster = arcpy.sa.CostDistance("SelectedPoint", SpeedRaster, None, None, None, None, None, None, '') 
%               out_distance_raster.save("CostPoint")
%               #out_raster = arcpy.sa.Divide("CostPoint", 1000)
%               #out_raster.save("Divide_Cost")
%               outfile = OutBase + str(row[0])+ ".asc"
%               arcpy.conversion.RasterToASCII("CostPoint", outfile)
%           #f = ["SelectedPoint","CostPoint","Divide_Cost"]
%           #for x in f:
%              # arcpy.Delete_management(x)