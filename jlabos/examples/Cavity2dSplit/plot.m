#!/usr/bin/octave
% Author : Yann Sagon <yann.sagon@unige.ch>

nx = 301;
ny = 301;

fp = load('res0.txt');
res = reshape(fp, nx, ny);
subplot (2, 2, 1)
imagesc(res);

fp = load('res1.txt');
res = reshape(fp, nx, ny);
subplot (2, 2, 2)
imagesc(res);

fp = load('res2.txt');
res = reshape(fp, nx, ny);
subplot (2, 2, 3)
imagesc(res);

fp = load('res3.txt');
res = reshape(fp, nx, ny);
subplot (2, 2, 4)
imagesc(res);

pause;
