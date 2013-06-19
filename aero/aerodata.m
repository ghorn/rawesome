function [ ret, data ] = aerodata
%AERODATA Summary of this function goes here
%   Detailed explanation goes here

% Alpha,pb/2V ,p'b/2V,Beta,qc/2V,Mach,rb/2V,r'b/2V,CX tot,Cl tot,Cl' tot,CY tot,Cm tot,CZ tot,Cn tot,Cn' tot,CL tot,CD tot,CD vis,CD ind,CL ff,CD ff,CY ff,e,Aileron,Elevator,CL a,CL b,CY a,CY b,Cl a,Cl b,Cm a,Cm b,Cn a,Cn b,CL p,CL q,CL r,CY p,CY q,CY r,Cl p,Cl q,Cl r,Cm p,Cm q,Cm r,Cn p,Cn q,Cn r,CL d1,CL d2,CY d1,CY d2,Cl d1,Cl d2,Cm d1,Cm d2,Cn d1,Cn d2,CD ff d1,CD ff d2,e d1,e d2,Xnp,Clb Cnr / Clr Cnb


% load('betty_alpFlaps.mat'); data=avltemphack_cdflaps;
% load('betty_alpElev.mat');  data=avltemphack_cdelev;
% load('betty_betAil.mat');   data=avltemphack_cdail;
% load('betty_betRud.mat');   data=avltemphack_cdrud;
% load('data_betty_jun7.mat'); data=avltemphack;
 load('ariane_ext.mat'); data=ariane_ext;


ret = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uncomment for Betty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ret.Alpha              = data(:,01);
% %ret.pb/2V             = data(:,02);
% %ret.p'b/2V            = data(:,03);
% ret.Beta               = data(:,04);
% %ret.qc/2V             = data(:,05);
% ret.Mach               = data(:,06);
% %ret.rb/2V             = data(:,07);
% %ret.r'b/2V            = data(:,08);
% ret.CX_tot             = data(:,09);
% 
% ret.Cl_tot             = data(:,10);
% %ret.Cl'_tot           = data(:,11);
% ret.CY_tot             = data(:,12);
% ret.Cm_tot             = data(:,13);
% ret.CZ_tot             = data(:,14);
% ret.Cn_tot             = data(:,15);
% %ret.Cn'_tot           = data(:,16);
% ret.CL_tot             = data(:,17);
% ret.CD_tot             = data(:,18);
% ret.CD_vis             = data(:,19);
% 
% ret.CD_ind             = data(:,20);
% ret.CL_ff              = data(:,21);
% ret.CD_ff              = data(:,22);
% ret.CY_ff              = data(:,23);
% ret.e                  = data(:,24);
% ret.Flap               = data(:,25);
% ret.Aileron            = data(:,26);
% ret.Elevator           = data(:,27);
% ret.Rudder             = data(:,28);
% 
% ret.CL_a               = data(:,29);
% ret.CL_b               = data(:,30);
% ret.CY_a               = data(:,31);
% 
% ret.CY_b               = data(:,32);
% ret.Cl_a               = data(:,33);
% ret.Cl_b               = data(:,34);
% ret.Cm_a               = data(:,35);
% ret.Cm_b               = data(:,36);
% ret.Cn_a               = data(:,37);
% ret.Cn_b               = data(:,38);
% ret.CL_p               = data(:,39);
% ret.CL_q               = data(:,40);
% ret.CL_r               = data(:,41);
% 
% ret.CY_p               = data(:,42);
% ret.CY_q               = data(:,43);
% ret.CY_r               = data(:,44);
% ret.Cl_p               = data(:,45);
% ret.Cl_q               = data(:,46);
% ret.Cl_r               = data(:,47);
% ret.Cm_p               = data(:,48);
% ret.Cm_q               = data(:,49);
% ret.Cm_r               = data(:,50);
% ret.Cn_p               = data(:,51);
% ret.Cn_q               = data(:,52);
% ret.Cn_r               = data(:,53);
% 
% ret.CL_d1              = data(:,54);
% ret.CL_d2              = data(:,55);
% ret.CL_d3              = data(:,56);
% ret.CL_d4              = data(:,57);
% 
% ret.CY_d1              = data(:,58);
% ret.CY_d2              = data(:,59);
% ret.CY_d3              = data(:,60);
% ret.CY_d4              = data(:,61);
% 
% ret.Cl_d1              = data(:,62);
% ret.Cl_d2              = data(:,63);
% ret.Cl_d3              = data(:,64);
% ret.Cl_d4              = data(:,65);
% 
% ret.Cm_d1              = data(:,66);
% ret.Cm_d2              = data(:,67);
% ret.Cm_d3              = data(:,68);
% ret.Cm_d4              = data(:,69);
% 
% ret.Cn_d1              = data(:,70);
% ret.Cn_d2              = data(:,71);
% ret.Cn_d3              = data(:,72);
% ret.Cn_d4              = data(:,73);
% 
% ret.CD_ff_d1           = data(:,74);
% ret.CD_ff_d2           = data(:,75);
% ret.CD_ff_d3           = data(:,76);
% ret.CD_ff_d4           = data(:,77);
% 
% ret.e_d1               = data(:,78);
% ret.e_d2               = data(:,79);
% ret.e_d3               = data(:,80);
% ret.e_d4               = data(:,81);
% 
% ret.Xnp                = data(:,82);
%ret.Clb_Cnr / Clr_Cnb  = data(:,75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uncomment for Ariane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ret.Alpha              = data(:,01);
%ret.pb/2V             = data(:,02);
%ret.p'b/2V            = data(:,03);
ret.Beta               = data(:,04);
%ret.qc/2V             = data(:,05);
ret.Mach               = data(:,06);
%ret.rb/2V             = data(:,07);
%ret.r'b/2V            = data(:,08);
ret.CX_tot             = data(:,09);

ret.Cl_tot             = data(:,10);
%ret.Cl'_tot           = data(:,11);
ret.CY_tot             = data(:,12);
ret.Cm_tot             = data(:,13);
ret.CZ_tot             = data(:,14);
ret.Cn_tot             = data(:,15);
%ret.Cn'_tot           = data(:,16);
ret.CL_tot             = data(:,17);
ret.CD_tot             = data(:,18);
ret.CD_vis             = data(:,19);

ret.CD_ind             = data(:,20);
ret.CL_ff              = data(:,21);
ret.CD_ff              = data(:,22);
ret.CY_ff              = data(:,23);
ret.e                  = data(:,24);
ret.Aileron            = data(:,25);
ret.Elevator           = data(:,26);

ret.CL_a               = data(:,27);
ret.CL_b               = data(:,28);
ret.CY_a               = data(:,29);

ret.CY_b               = data(:,30);
ret.Cl_a               = data(:,31);
ret.Cl_b               = data(:,32);
ret.Cm_a               = data(:,33);
ret.Cm_b               = data(:,34);
ret.Cn_a               = data(:,35);
ret.Cn_b               = data(:,36);
ret.CL_p               = data(:,37);
ret.CL_q               = data(:,38);
ret.CL_r               = data(:,39);

ret.CY_p               = data(:,40);
ret.CY_q               = data(:,41);
ret.CY_r               = data(:,42);
ret.Cl_p               = data(:,43);
ret.Cl_q               = data(:,44);
ret.Cl_r               = data(:,45);
ret.Cm_p               = data(:,46);
ret.Cm_q               = data(:,47);
ret.Cm_r               = data(:,48);
ret.Cn_p               = data(:,49);
ret.Cn_q               = data(:,50);
ret.Cn_r               = data(:,51);

ret.CL_d1              = data(:,52);
ret.CL_d2              = data(:,53);

ret.CY_d1              = data(:,54);
ret.CY_d2              = data(:,55);

ret.Cl_d1              = data(:,56);
ret.Cl_d2              = data(:,57);

ret.Cm_d1              = data(:,58);
ret.Cm_d2              = data(:,59);

ret.Cn_d1              = data(:,60);
ret.Cn_d2              = data(:,61);

ret.CD_ff_d1           = data(:,62);
ret.CD_ff_d2           = data(:,63);

ret.e_d1               = data(:,64);
ret.e_d2               = data(:,65);

ret.Xnp                = data(:,66);
end
