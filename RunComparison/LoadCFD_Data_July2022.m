function [Lam,DNS,RANS] = LoadCFD_Data_July2022(cfd_filepath,Run_Mean_Velo)
    
Lam     = struct;
DNS     = struct;
RANS    = struct;

    DNS.Velo = readtable(cfd_filepath,'Range','E2:E202');
    DNS.height = readtable(cfd_filepath,'Range','J2:J202');
    Lam.Velo = readtable(cfd_filepath,'Range','O3:O203');
    Lam.height = readtable(cfd_filepath,'Range','Q3:Q203');
    RANS.Velo = readtable(cfd_filepath,'Range','V3:V203');
    RANS.height = readtable(cfd_filepath,'Range','X3:X203');
    DNS.Velo = DNS.Velo{:,:};
    DNS.height = DNS.height{:,:};
    Lam.Velo = Lam.Velo{:,:};
    Lam.height = Lam.height{:,:};
    RANS.Velo = RANS.Velo{:,:};
    RANS.height = RANS.height{:,:};
    scale_velo = mean(Run_Mean_Velo(:,1))./max(Lam.Velo);
    DNS.Velo = scale_velo.*DNS.Velo;
    Lam.Velo = scale_velo.*Lam.Velo;
    RANS.Velo = scale_velo.*RANS.Velo;

end