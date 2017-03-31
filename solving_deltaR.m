%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Martin Rodriguez Cruz
% Solving delta R
% Created 2/21/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc


syms a b s00 s01 s10 s11;
assume(a,'real'); assume(b,'real');
assume(s00,'real'); assume(s01,'real');
assume(s10,'real'); assume(s11,'real');

RTdR = [ 0 a; -a 0 ] 
LHS = [ 0 b; -b 0 ]
S = [s00, s01; s10, s11]

A = RTdR*S + S*RTdR
A(1,2)

[M, B] = equationsToMatrix( [A(1,2) == b],[a] ) 

% syms a b c d e f s00 s01 s02 s10 s11 s12 s20 s21 s22
% assume(a,'real'); assume(b,'real'); assume(c,'real');
% assume(d,'real'); assume(e,'real'); assume(f,'real');
% assume(s00,'real'); assume(s01,'real'); assume(s02,'real');
% assume(s10,'real'); assume(s11,'real'); assume(s12,'real');
% assume(s20,'real'); assume(s21,'real'); assume(s22,'real');
% 
% 
% RTdR = [0 a b; -a 0 c; -b -c 0];
% LHS = [0 d e; -d 0 f; -e -f 0];
% S = [s00 s01 s02; s10 s11 s12; s20 s21 s22];
% 
% A = RTdR*S + S*RTdR
% A(1,2)
% A(1,3)
% A(2,3)
% 
% [M, b] = equationsToMatrix( [A(1,2) == d, A(1,3) == e, A(2,3) == f],[a, b, c] )