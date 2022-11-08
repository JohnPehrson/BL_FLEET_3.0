function [Lam,PB,SRA] = Load_July2022_Data(dist_LE_Trips,FLEET_July_2022_filepath)
    Lam = struct;
    Lam.Mean_Velo       = readmatrix(FLEET_July_2022_filepath,'Sheet','Laminar Array','Range','A2:A323'); %seconds
    Lam.Mean_Velo_Unc   = readmatrix(FLEET_July_2022_filepath,'Sheet','Laminar Array','Range','B2:B323'); %seconds
    Lam.Heights         = readmatrix(FLEET_July_2022_filepath,'Sheet','Laminar Array','Range','J2:J323'); %seconds
    Lam.Downstream_Loc  = dist_LE_Trips + readmatrix(FLEET_July_2022_filepath,'Sheet','Laminar Array','Range','F2:F2'); %seconds
    
    PB = struct;
    PB.Mean_Velo       = readmatrix(FLEET_July_2022_filepath,'Sheet','Pizza Box Array','Range','A2:A323'); %seconds
    PB.Mean_Velo_Unc   = readmatrix(FLEET_July_2022_filepath,'Sheet','Pizza Box Array','Range','B2:B323'); %seconds
    PB.Heights         = readmatrix(FLEET_July_2022_filepath,'Sheet','Pizza Box Array','Range','J2:J323'); %seconds
    PB.Downstream_Loc  = dist_LE_Trips + readmatrix(FLEET_July_2022_filepath,'Sheet','Pizza Box Array','Range','F2:F2'); %seconds
    
    SRA = struct;
    SRA.Mean_Velo       = readmatrix(FLEET_July_2022_filepath,'Sheet','Synthetic Roughness Array','Range','A2:A323'); %seconds
    SRA.Mean_Velo_Unc   = readmatrix(FLEET_July_2022_filepath,'Sheet','Synthetic Roughness Array','Range','B2:B323'); %seconds
    SRA.Heights         = readmatrix(FLEET_July_2022_filepath,'Sheet','Synthetic Roughness Array','Range','J2:J323'); %seconds
    SRA.Downstream_Loc  = dist_LE_Trips + readmatrix(FLEET_July_2022_filepath,'Sheet','Synthetic Roughness Array','Range','F2:F2'); %seconds
end