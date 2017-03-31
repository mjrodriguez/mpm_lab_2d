%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Martin Rodriguez Cruz
% Data Visualizer for Particles
% Created 2/14/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


syms F00 F01 F10 F11;
assume(F00,'real');
assume(F01,'real');
assume(F10,'real');
assume(F11,'real');

F = [ F00, F01; F10 F11 ]
J = det(F)
dJFT = cell(2,2); 
for i = 1:2
    for j = 1:2
       
        disp(i)
        disp(j)
        
       dJFT{i,j} = diff(J*inv(F)',F(i,j));
       disp(dJFT{i,j});
       
    end
end

% syms F00 F01 F02 F10 F11 F12 F20 F21 F22;
% assume(F00,'real');
% assume(F01,'real');
% assume(F02,'real');
% assume(F10,'real');
% assume(F11,'real');
% assume(F12,'real');
% assume(F20,'real');
% assume(F21,'real');
% assume(F22,'real');
% 
% F = [F00 F01 F02; F10 F11 F12; F20 F21 F22];
% J = det(F);
% 
% dJFT = cell(3,3); 
% for i = 1:3
%     for j = 1:3
%        
%         disp(i)
%         disp(j)
%         
%        dJFT{i,j} = diff(J*inv(F)',F(i,j));
%        disp(dJFT{i,j});
%        
%     end
% end




