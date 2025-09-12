{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica-Bold;\f1\fswiss\fcharset0 Helvetica;\f2\fswiss\fcharset0 Helvetica-Oblique;
}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\b\fs24 \cf0 Some instructions about the code and the variables here\

\f1\b0 \
I will give some examples of plotting code for the given data file, but this will be commented out, and you can uncomment it as needed.\
\

\f2\i main_fm 
\f1\i0 is the main function script and here is everything you need to run the code and the variables that you can change, extract etc.\
\
Some variables in 
\f2\i main_fm 
\f1\i0  and what their function is:\
\
Input - 
\f2\i match_travletimes_mean_trend 
\f1\i0  \
*! First, in order to match the traveltimes with corresponding prior mean, input traveltimes in the function. Pass it as a matrix. [Let me know if other functionality us preferred).\
\
Output - 
\f2\i match_travletimes_mean_trend
\f1\i0 \
- 
\f2\i prior_mean: 
\f1\i0  dim: ((n x m) x 4); col 1 - x_g, col 2 - x_o, col 3 - x_cl, col 4 - depth\
\
Input: 
\f2\i sim_forward_model
\f1\i0 \
\
- 
\f2\i x_g: 
\f1\i0 a vector of gas saturations - saturation, Gaussian\
- 
\f2\i x_o: 
\f1\i0 a vector of oil saturations - saturation, Gaussian\
- 
\f2\i x_cl : 
\f1\i0 a vector of clay content - Gaussian \
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 - 
\f2\i traveltimes: 
\f1\i0 a section of the traveltimes matrix passed as argument, either a single location or a matrix, or a vector (whatever is passed, is transformed to appropriate shape\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 \
Shape of the output will be:\
\
- 
\f2\i  vector of dimension ((n x m)x2), 
\f1\i0  format being (R_0,G)_1, (R_0,G)_2,\'85. \
\
\
For plotting of traveltimes/depth:\
\
Extract the preferred area from the traveltimes matrix AND extract the corresponding area from inline and from xline matrix. Plot this using for example image.plot. \
\
\
\
\
}
