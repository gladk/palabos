#!/usr/bin/octave
% Author : Yann Sagon <yann.sagon@unige.ch>

nx = 121;
ny = 801;

fp = load('res0.txt');
res = reshape(fp, nx, ny);
subplot (3, 1, 1)
imagesc(res);

fp = load('res1.txt');
res = reshape(fp, nx, ny);
subplot (3, 1, 2)
imagesc(res);

fp = load('res2.txt');
res = reshape(fp, nx, ny);
subplot (3, 1, 3)
imagesc(res);

pause;
