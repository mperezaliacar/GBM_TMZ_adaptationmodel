clear all
close all
clc

%% cálculo del pvalor

chi2sim = 7.919;
nu = 22;

p = 1-chi2cdf(chi2sim,nu)