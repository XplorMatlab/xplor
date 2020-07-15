function c = vdaqcolors
% function c = vdaqcolors
%---
% returns VDAQ color map

c = [hex2dec(['FF';'FF';'FF'])'
   hex2dec(['FF';'00';'00'])'
   hex2dec(['FF';'00';'CC'])'
   hex2dec(['66';'00';'66'])'
   hex2dec(['66';'66';'99'])'
   hex2dec(['00';'00';'ff'])'
   hex2dec(['00';'00';'66'])'
   hex2dec(['FF';'FF';'00'])'
   hex2dec(['00';'FF';'00'])'
   hex2dec(['66';'60';'00'])'
   hex2dec(['00';'33';'33'])'
   hex2dec(['FF';'99';'00'])'
   hex2dec(['FF';'33';'00'])'
   hex2dec(['cc';'00';'00'])'
   hex2dec(['66';'33';'00'])'];
c = c(end:-1:1,:)/255; 