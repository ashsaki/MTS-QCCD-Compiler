OPENQASM 2.0;
include "qelib1.inc";
qreg q1[64];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
cx q1[47],q1[54];
cx q1[47],q1[40];
cx q1[47],q1[38];
cx q1[54],q1[33];
cx q1[54],q1[50];
cx q1[39],q1[59];
cx q1[39],q1[8];
cx q1[39],q1[13];
cx q1[59],q1[24];
cx q1[59],q1[41];
cx q1[7],q1[25];
cx q1[7],q1[18];
cx q1[7],q1[45];
cx q1[25],q1[16];
cx q1[25],q1[3];
cx q1[37],q1[48];
cx q1[37],q1[44];
cx q1[37],q1[5];
cx q1[48],q1[36];
cx q1[48],q1[51];
cx q1[9],q1[21];
cx q1[9],q1[58];
cx q1[9],q1[53];
cx q1[21],q1[30];
cx q1[21],q1[60];
cx q1[44],q1[56];
cx q1[44],q1[14];
cx q1[56],q1[38];
cx q1[56],q1[11];
cx q1[16],q1[61];
cx q1[16],q1[36];
cx q1[2],q1[6];
cx q1[2],q1[61];
cx q1[2],q1[58];
cx q1[6],q1[23];
cx q1[6],q1[22];
cx q1[13],q1[35];
cx q1[13],q1[14];
cx q1[35],q1[49];
cx q1[35],q1[63];
cx q1[29],q1[30];
cx q1[29],q1[10];
cx q1[29],q1[57];
cx q1[30],q1[10];
cx q1[23],q1[53];
cx q1[23],q1[28];
cx q1[53],q1[12];
cx q1[38],q1[1];
cx q1[28],q1[49];
cx q1[28],q1[43];
cx q1[11],q1[41];
cx q1[11],q1[32];
cx q1[41],q1[5];
cx q1[36],q1[55];
cx q1[55],q1[52];
cx q1[55],q1[46];
cx q1[1],q1[62];
cx q1[1],q1[27];
cx q1[10],q1[12];
cx q1[5],q1[42];
cx q1[24],q1[17];
cx q1[24],q1[45];
cx q1[4],q1[31];
cx q1[4],q1[18];
cx q1[4],q1[34];
cx q1[31],q1[20];
cx q1[31],q1[32];
cx q1[33],q1[40];
cx q1[33],q1[17];
cx q1[40],q1[52];
cx q1[17],q1[19];
cx q1[49],q1[3];
cx q1[18],q1[57];
cx q1[57],q1[52];
cx q1[20],q1[22];
cx q1[20],q1[34];
cx q1[0],q1[8];
cx q1[0],q1[50];
cx q1[0],q1[26];
cx q1[8],q1[27];
cx q1[46],q1[63];
cx q1[46],q1[19];
cx q1[63],q1[45];
cx q1[22],q1[43];
cx q1[43],q1[32];
cx q1[3],q1[15];
cx q1[15],q1[42];
cx q1[15],q1[60];
cx q1[26],q1[61];
cx q1[26],q1[19];
cx q1[51],q1[58];
cx q1[51],q1[14];
cx q1[50],q1[42];
cx q1[60],q1[12];
cx q1[62],q1[27];
cx q1[62],q1[34];
