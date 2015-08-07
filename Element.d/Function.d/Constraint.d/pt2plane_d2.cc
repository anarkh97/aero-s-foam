#include <math.h>

void pt2plane_d2 (double x[12], double ddd_dx_dx[144])
{
  double d[1];
  double t1;
  double t10;
  double t100;
  double t105;
  double t11;
  double t110;
  double t113;
  double t117;
  double t12;
  double t13;
  double t14;
  double t142;
  double t15;
  double t16;
  double t169;
  double t171;
  double t174;
  double t175;
  double t178;
  double t179;
  double t18;
  double t180;
  double t188;
  double t19;
  double t2;
  double t20;
  double t201;
  double t21;
  double t214;
  double t22;
  double t228;
  double t23;
  double t233;
  double t24;
  double t246;
  double t25;
  double t251;
  double t26;
  double t267;
  double t28;
  double t281;
  double t286;
  double t29;
  double t299;
  double t3;
  double t30;
  double t304;
  double t31;
  double t318;
  double t32;
  double t321;
  double t322;
  double t34;
  double t35;
  double t36;
  double t363;
  double t37;
  double t38;
  double t39;
  double t4;
  double t40;
  double t404;
  double t419;
  double t42;
  double t422;
  double t423;
  double t43;
  double t44;
  double t46;
  double t49;
  double t5;
  double t511;
  double t514;
  double t515;
  double t52;
  double t55;
  double t579;
  double t58;
  double t582;
  double t583;
  double t6;
  double t61;
  double t638;
  double t64;
  double t641;
  double t642;
  double t67;
  double t688;
  double t691;
  double t692;
  double t7;
  double t70;
  double t716;
  double t719;
  double t720;
  double t735;
  double t738;
  double t739;
  double t74;
  double t79;
  double t8;
  double t82;
  double t85;
  double t9;
  double t90;
  double t94;
  double t97;
  t1 = x[7];
  t2 = x[4];
  t3 = t1 - t2;
  t4 = x[11];
  t5 = x[5];
  t6 = t4 - t5;
  t7 = t3 * t6;
  t8 = x[8];
  t9 = t8 - t5;
  t10 = x[10];
  t11 = t10 - t2;
  t12 = t9 * t11;
  t13 = t7 - t12;
  t14 = x[0];
  t15 = x[3];
  t16 = t14 - t15;
  t18 = x[6];
  t19 = t18 - t15;
  t20 = t19 * t6;
  t21 = x[9];
  t22 = t21 - t15;
  t23 = t9 * t22;
  t24 = -t20 + t23;
  t25 = x[1];
  t26 = t25 - t2;
  t28 = t19 * t11;
  t29 = t3 * t22;
  t30 = t28 - t29;
  t31 = x[2];
  t32 = t31 - t5;
  t34 = t13 * t16 + t24 * t26 + t30 * t32;
  t35 = t13 * t13;
  t36 = t24 * t24;
  t37 = t30 * t30;
  t38 = t35 + t36 + t37;
  t39 = sqrt(t38);
  t40 = 1.0 / t39;
  d[0] = t34 * t40;
  ddd_dx_dx[0] = 0.0;
  ddd_dx_dx[1] = 0.0;
  ddd_dx_dx[2] = 0.0;
  t42 = 1.0 / t39 / t38;
  t43 = t13 * t42;
  t44 = t4 - t8;
  t46 = -t10 + t1;
  t49 = 2.0 * t24 * t44 + 2.0 * t30 * t46;
  ddd_dx_dx[3] = -t43 * t49 / 2.0;
  t52 = -t44;
  t55 = -t18 + t21;
  t58 = 2.0 * t13 * t52 + 2.0 * t30 * t55;
  ddd_dx_dx[4] = t52 * t40 - t43 * t58 / 2.0;
  t61 = -t46;
  t64 = -t55;
  t67 = 2.0 * t13 * t61 + 2.0 * t24 * t64;
  ddd_dx_dx[5] = t61 * t40 - t43 * t67 / 2.0;
  t70 = -t6;
  t74 = 2.0 * t11 * t30 + 2.0 * t24 * t70;
  ddd_dx_dx[6] = -t43 * t74 / 2.0;
  t79 = -t22;
  t82 = 2.0 * t13 * t6 + 2.0 * t30 * t79;
  ddd_dx_dx[7] = t6 * t40 - t43 * t82 / 2.0;
  t85 = -t11;
  t90 = 2.0 * t13 * t85 + 2.0 * t22 * t24;
  ddd_dx_dx[8] = t85 * t40 - t43 * t90 / 2.0;
  t94 = -t3;
  t97 = 2.0 * t24 * t9 + 2.0 * t30 * t94;
  ddd_dx_dx[9] = -t43 * t97 / 2.0;
  t100 = -t9;
  t105 = 2.0 * t100 * t13 + 2.0 * t19 * t30;
  ddd_dx_dx[10] = t100 * t40 - t43 * t105 / 2.0;
  t110 = -t19;
  t113 = 2.0 * t110 * t24 + 2.0 * t13 * t3;
  ddd_dx_dx[11] = t3 * t40 - t43 * t113 / 2.0;
  ddd_dx_dx[12] = 0.0e0;
  ddd_dx_dx[13] = 0.0e0;
  ddd_dx_dx[14] = 0.0e0;
  t117 = t24 * t42;
  ddd_dx_dx[15] = t44 * t40 - t117 * t49 / 2.0;
  ddd_dx_dx[16] = -t117 * t58 / 2.0;
  ddd_dx_dx[17] = t64 * t40 - t117 * t67 / 2.0;
  ddd_dx_dx[18] = t70 * t40 - t117 * t74 / 2.0;
  ddd_dx_dx[19] = -t117 * t82 / 2.0;
  ddd_dx_dx[20] = t22 * t40 - t117 * t90 / 2.0;
  ddd_dx_dx[21] = t9 * t40 - t117 * t97 / 2.0;
  ddd_dx_dx[22] = -t117 * t105 / 2.0;
  ddd_dx_dx[23] = t110 * t40 - t117 * t113 / 2.0;
  ddd_dx_dx[24] = 0.0e0;
  ddd_dx_dx[25] = 0.0e0;
  ddd_dx_dx[26] = 0.0e0;
  t142 = t30 * t42;
  ddd_dx_dx[27] = t46 * t40 - t142 * t49 / 2.0;
  ddd_dx_dx[28] = t55 * t40 - t142 * t58 / 2.0;
  ddd_dx_dx[29] = -t142 * t67 / 2.0;
  ddd_dx_dx[30] = t11 * t40 - t142 * t74 / 2.0;
  ddd_dx_dx[31] = t79 * t40 - t142 * t82 / 2.0;
  ddd_dx_dx[32] = -t142 * t90 / 2.0;
  ddd_dx_dx[33] = t94 * t40 - t142 * t97 / 2.0;
  ddd_dx_dx[34] = t19 * t40 - t142 * t105 / 2.0;
  ddd_dx_dx[35] = -t142 * t113 / 2.0;
  ddd_dx_dx[36] = ddd_dx_dx[3];
  ddd_dx_dx[37] = ddd_dx_dx[15];
  ddd_dx_dx[38] = ddd_dx_dx[27];
  t169 = (t26 * t44 + t32 * t46 + t12 - t7) * t42;
  t171 = t38 * t38;
  t174 = t34 / t39 / t171;
  t175 = t49 * t49;
  t178 = t34 * t42;
  t179 = t44 * t44;
  t180 = t46 * t46;
  ddd_dx_dx[39] = -t169 * t49 + 3.0 / 4.0 * t174 * t175 - t178 * (2.0 * t179 + 2.0 * t180) / 2.0;
  t188 = (t16 * t52 + t32 * t55 + t20 - t23) * t42;
  ddd_dx_dx[40] = -t188 * t49 / 2.0 - t169 * t58 / 2.0 + 3.0 / 4.0 * t174 * t58 * t49 - t178 * t46 * t55;
  t201 = (t16 * t61 + t26 * t64 - t28 + t29) * t42;
  ddd_dx_dx[41] = -t201 * t49 / 2.0 - t169 * t67 / 2.0 + 3.0 / 4.0 * t174 * t67 * t49 - t178 * t44 * t64;
  t214 = (t11 * t32 + t26 * t70) * t42;
  ddd_dx_dx[42] = -t214 * t49 / 2.0 - t169 * t74 / 2.0 + 3.0 / 4.0 * t174 * t74 * t49 - t178 * (2.0 * t11 * t46 + 2.0 * t44 * t70) / 2.0;
  t228 = -t4 + t31;
  t233 = (t16 * t6 + t32 * t79) * t42;
  ddd_dx_dx[43] = t228 * t40 - t233 * t49 / 2.0 - t169 * t82 / 2.0 + 3.0 / 4.0 * t174 * t82 * t49 - t178 * (2.0 * t46 * t79 + 2.0 * t28 - 2.0 * t29) / 2.0;
  t246 = t10 - t25;
  t251 = (t16 * t85 + t22 * t26) * t42;
  ddd_dx_dx[44] = t246 * t40 - t251 * t49 / 2.0 - t169 * t90 / 2.0 + 3.0 / 4.0 * t174 * t90 * t49 - t178 * (2.0 * t22 * t44 + 2.0 * t20 - 2.0 * t23) / 2.0;
  t267 = (t26 * t9 + t32 * t94) * t42;
  ddd_dx_dx[45] = -t267 * t49 / 2.0 - t169 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t49 - t178 * (2.0 * t44 * t9 + 2.0 * t46 * t94) / 2.0;
  t281 = t8 - t31;
  t286 = (t100 * t16 + t19 * t32) * t42;
  ddd_dx_dx[46] = t281 * t40 - t286 * t49 / 2.0 - t169 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t49 - t178 * (2.0 * t19 * t46 - 2.0 * t28 + 2.0 * t29) / 2.0;
  t299 = -t1 + t25;
  t304 = (t110 * t26 + t16 * t3) * t42;
  ddd_dx_dx[47] = t299 * t40 - t304 * t49 / 2.0 - t169 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t49 - t178 * (2.0 * t110 * t44 - 2.0 * t20 + 2.0 * t23) / 2.0;
  ddd_dx_dx[48] = ddd_dx_dx[4];
  ddd_dx_dx[49] = ddd_dx_dx[16];
  ddd_dx_dx[50] = ddd_dx_dx[28];
  ddd_dx_dx[51] = ddd_dx_dx[40];
  t318 = t58 * t58;
  t321 = t52 * t52;
  t322 = t55 * t55;
  ddd_dx_dx[52] = -t188 * t58 + 3.0 / 4.0 * t174 * t318 - t178 * (2.0 * t321 + 2.0 * t322) / 2.0;
  ddd_dx_dx[53] = -t201 * t58 / 2.0 - t188 * t67 / 2.0 + 3.0 / 4.0 * t174 * t67 * t58 - t178 * t52 * t61;
  ddd_dx_dx[54] = -t228 * t40 - t214 * t58 / 2.0 - t188 * t74 / 2.0 + 3.0 / 4.0 * t174 * t74 * t58 - t178 * (2.0 * t11 * t55 - 2.0 * t28 + 2.0 * t29) / 2.0;
  ddd_dx_dx[55] = -t233 * t58 / 2.0 - t188 * t82 / 2.0 + 3.0 / 4.0 * t174 * t82 * t58 - t178 * (2.0 * t52 * t6 + 2.0 * t55 * t79) / 2.0;
  t363 = t14 - t21;
  ddd_dx_dx[56] = t363 * t40 - t251 * t58 / 2.0 - t188 * t90 / 2.0 + 3.0 / 4.0 * t174 * t90 * t58 - t178 * (2.0 * t52 * t85 - 2.0 * t12 + 2.0 * t7) / 2.0;
  ddd_dx_dx[57] = -t281 * t40 - t267 * t58 / 2.0 - t188 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t58 - t178 * (2.0 * t55 * t94 + 2.0 * t28 - 2.0 * t29) / 2.0;
  ddd_dx_dx[58] = -t286 * t58 / 2.0 - t188 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t58 - t178 * (2.0 * t100 * t52 + 2.0 * t19 * t55) / 2.0;
  t404 = -t14 + t18;
  ddd_dx_dx[59] = t404 * t40 - t304 * t58 / 2.0 - t188 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t58 - t178 * (2.0 * t3 * t52 + 2.0 * t12 - 2.0 * t7) / 2.0;
  ddd_dx_dx[60] = ddd_dx_dx[5];
  ddd_dx_dx[61] = ddd_dx_dx[17];
  ddd_dx_dx[62] = ddd_dx_dx[29];
  ddd_dx_dx[63] = ddd_dx_dx[41];
  ddd_dx_dx[64] = ddd_dx_dx[53];
  t419 = t67 * t67;
  t422 = t61 * t61;
  t423 = t64 * t64;
  ddd_dx_dx[65] = -t201 * t67 + 3.0 / 4.0 * t174 * t419 - t178 * (2.0 * t422 + 2.0 * t423) / 2.0;
  ddd_dx_dx[66] = -t246 * t40 - t214 * t67 / 2.0 - t201 * t74 / 2.0 + 3.0 / 4.0 * t174 * t74 * t67 - t178 * (2.0 * t64 * t70 - 2.0 * t20 + 2.0 * t23) / 2.0;
  ddd_dx_dx[67] = -t363 * t40 - t233 * t67 / 2.0 - t201 * t82 / 2.0 + 3.0 / 4.0 * t174 * t82 * t67 - t178 * (2.0 * t6 * t61 + 2.0 * t12 - 2.0 * t7) / 2.0;
  ddd_dx_dx[68] = -t251 * t67 / 2.0 - t201 * t90 / 2.0 + 3.0 / 4.0 * t174 * t90 * t67 - t178 * (2.0 * t22 * t64 + 2.0 * t61 * t85) / 2.0;
  ddd_dx_dx[69] = -t299 * t40 - t267 * t67 / 2.0 - t201 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t67 - t178 * (2.0 * t64 * t9 + 2.0 * t20 - 2.0 * t23) / 2.0;
  ddd_dx_dx[70] = -t404 * t40 - t286 * t67 / 2.0 - t201 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t67 - t178 * (2.0 * t100 * t61 - 2.0 * t12 + 2.0 * t7) / 2.0;
  ddd_dx_dx[71] = -t304 * t67 / 2.0 - t201 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t67 - t178 * (2.0 * t110 * t64 + 2.0 * t3 * t61) / 2.0;
  ddd_dx_dx[72] = ddd_dx_dx[6];
  ddd_dx_dx[73] = ddd_dx_dx[18];
  ddd_dx_dx[74] = ddd_dx_dx[30];
  ddd_dx_dx[75] = ddd_dx_dx[42];
  ddd_dx_dx[76] = ddd_dx_dx[54];
  ddd_dx_dx[77] = ddd_dx_dx[66];
  t511 = t74 * t74;
  t514 = t70 * t70;
  t515 = t11 * t11;
  ddd_dx_dx[78] = -t214 * t74 + 3.0 / 4.0 * t174 * t511 - t178 * (2.0 * t514 + 2.0 * t515) / 2.0;
  ddd_dx_dx[79] = -t233 * t74 / 2.0 - t214 * t82 / 2.0 + 3.0 / 4.0 * t174 * t82 * t74 - t178 * t11 * t79;
  ddd_dx_dx[80] = -t251 * t74 / 2.0 - t214 * t90 / 2.0 + 3.0 / 4.0 * t174 * t90 * t74 - t178 * t70 * t22;
  ddd_dx_dx[81] = -t267 * t74 / 2.0 - t214 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t74 - t178 * (2.0 * t11 * t94 + 2.0 * t70 * t9) / 2.0;
  ddd_dx_dx[82] = t32 * t40 - t286 * t74 / 2.0 - t214 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t74 - t178 * (4.0 * t28 - 2.0 * t29) / 2.0;
  ddd_dx_dx[83] = -t26 * t40 - t304 * t74 / 2.0 - t214 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t74 - t178 * (2.0 * t110 * t70 + 2.0 * t20 - 2.0 * t23) / 2.0;
  ddd_dx_dx[84] = ddd_dx_dx[7];
  ddd_dx_dx[85] = ddd_dx_dx[19];
  ddd_dx_dx[86] = ddd_dx_dx[31];
  ddd_dx_dx[87] = ddd_dx_dx[43];
  ddd_dx_dx[88] = ddd_dx_dx[55];
  ddd_dx_dx[89] = ddd_dx_dx[67];
  ddd_dx_dx[90] = ddd_dx_dx[79];
  t579 = t82 * t82;
  t582 = t6 * t6;
  t583 = t79 * t79;
  ddd_dx_dx[91] = -t233 * t82 + 3.0 / 4.0 * t174 * t579 - t178 * (2.0 * t582 + 2.0 * t583) / 2.0;
  ddd_dx_dx[92] = -t251 * t82 / 2.0 - t233 * t90 / 2.0 + 3.0 / 4.0 * t174 * t90 * t82 - t178 * t6 * t85;
  ddd_dx_dx[93] = -t32 * t40 - t267 * t82 / 2.0 - t233 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t82 - t178 * (2.0 * t79 * t94 - 2.0 * t28 + 2.0 * t29) / 2.0;
  ddd_dx_dx[94] = -t286 * t82 / 2.0 - t233 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t82 - t178 * (2.0 * t100 * t6 + 2.0 * t19 * t79) / 2.0;
  ddd_dx_dx[95] = t16 * t40 - t304 * t82 / 2.0 - t233 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t82 - t178 * (4.0 * t7 - 2.0 * t12) / 2.0;
  ddd_dx_dx[96] = ddd_dx_dx[8];
  ddd_dx_dx[97] = ddd_dx_dx[20];
  ddd_dx_dx[98] = ddd_dx_dx[32];
  ddd_dx_dx[99] = ddd_dx_dx[44];
  ddd_dx_dx[100] = ddd_dx_dx[56];
  ddd_dx_dx[101] = ddd_dx_dx[68];
  ddd_dx_dx[102] = ddd_dx_dx[80];
  ddd_dx_dx[103] = ddd_dx_dx[92];
  t638 = t90 * t90;
  t641 = t85 * t85;
  t642 = t22 * t22;
  ddd_dx_dx[104] = -t251 * t90 + 3.0 / 4.0 * t174 * t638 - t178 * (2.0 * t641 + 2.0 * t642) / 2.0;
  ddd_dx_dx[105] = t26 * t40 - t267 * t90 / 2.0 - t251 * t97 / 2.0 + 3.0 / 4.0 * t174 * t97 * t90 - t178 * (4.0 * t23 - 2.0 * t20) / 2.0;
  ddd_dx_dx[106] = -t16 * t40 - t286 * t90 / 2.0 - t251 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t90 - t178 * (2.0 * t100 * t85 + 2.0 * t12 - 2.0 * t7) / 2.0;
  ddd_dx_dx[107] = -t304 * t90 / 2.0 - t251 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t90 - t178 * (2.0 * t110 * t22 + 2.0 * t3 * t85) / 2.0;
  ddd_dx_dx[108] = ddd_dx_dx[9];
  ddd_dx_dx[109] = ddd_dx_dx[21];
  ddd_dx_dx[110] = ddd_dx_dx[33];
  ddd_dx_dx[111] = ddd_dx_dx[45];
  ddd_dx_dx[112] = ddd_dx_dx[57];
  ddd_dx_dx[113] = ddd_dx_dx[69];
  ddd_dx_dx[114] = ddd_dx_dx[81];
  ddd_dx_dx[115] = ddd_dx_dx[93];
  ddd_dx_dx[116] = ddd_dx_dx[105];
  t688 = t97 * t97;
  t691 = t9 * t9;
  t692 = t94 * t94;
  ddd_dx_dx[117] = -t267 * t97 + 3.0 / 4.0 * t174 * t688 - t178 * (2.0 * t691 + 2.0 * t692) / 2.0;
  ddd_dx_dx[118] = -t286 * t97 / 2.0 - t267 * t105 / 2.0 + 3.0 / 4.0 * t174 * t105 * t97 - t178 * t94 * t19;
  ddd_dx_dx[119] = -t304 * t97 / 2.0 - t267 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t97 - t178 * t9 * t110;
  ddd_dx_dx[120] = ddd_dx_dx[10];
  ddd_dx_dx[121] = ddd_dx_dx[22];
  ddd_dx_dx[122] = ddd_dx_dx[34];
  ddd_dx_dx[123] = ddd_dx_dx[46];
  ddd_dx_dx[124] = ddd_dx_dx[58];
  ddd_dx_dx[125] = ddd_dx_dx[70];
  ddd_dx_dx[126] = ddd_dx_dx[82];
  ddd_dx_dx[127] = ddd_dx_dx[94];
  ddd_dx_dx[128] = ddd_dx_dx[106];
  ddd_dx_dx[129] = ddd_dx_dx[118];
  t716 = t105 * t105;
  t719 = t100 * t100;
  t720 = t19 * t19;
  ddd_dx_dx[130] = -t286 * t105 + 3.0 / 4.0 * t174 * t716 - t178 * (2.0 * t719 + 2.0 * t720) / 2.0;
  ddd_dx_dx[131] = -t304 * t105 / 2.0 - t286 * t113 / 2.0 + 3.0 / 4.0 * t174 * t113 * t105 - t178 * t100 * t3;
  ddd_dx_dx[132] = ddd_dx_dx[11];
  ddd_dx_dx[133] = ddd_dx_dx[23];
  ddd_dx_dx[134] = ddd_dx_dx[35];
  ddd_dx_dx[135] = ddd_dx_dx[47];
  ddd_dx_dx[136] = ddd_dx_dx[59];
  ddd_dx_dx[137] = ddd_dx_dx[71];
  ddd_dx_dx[138] = ddd_dx_dx[83];
  ddd_dx_dx[139] = ddd_dx_dx[95];
  ddd_dx_dx[140] = ddd_dx_dx[107];
  ddd_dx_dx[141] = ddd_dx_dx[119];
  ddd_dx_dx[142] = ddd_dx_dx[131];
  t735 = t113 * t113;
  t738 = t3 * t3;
  t739 = t110 * t110;
  ddd_dx_dx[143] = -t304 * t113 + 3.0 / 4.0 * t174 * t735 - t178 * (2.0 * t738 + 2.0 * t739) / 2.0;
  return;
}
