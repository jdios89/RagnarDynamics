#include <Arduino.h>
#include "RagnarDynamics.h"
#include <Matrix.h>
#include <Ragnarkinematics.h>
#include <math.h>

void ragnarMassCentrGrav(
    double theta[4], double dtheta[4], float (*params)[4][8], 
    double parammass[6], double sc[4][6], double (*mass_matrix)[7][7], 
    double (*coriolis_vector)[7], double (*gravity_vector)[7], 
    double sct[8], double g)
{
    double sa1, ca1, sb1, cb1, sg1, cg1, sa2, ca2, sb2, cb2, sg2, cg2;
    double sa3, ca3, sb3, cb3, sg3, cg3, sa4, ca4, sb4, cb4, sg4, cg4;
    double st1, ct1, st2, ct2, st3, ct3, st4, ct4;
    double mp, mb, Ib, ml;
    // extract precomputed sines and cosines
    sa1 = sc[0][0];
    ca1 = sc[0][1];
    sb1 = sc[0][2];
    cb1 = sc[0][3];
    sg1 = sc[0][4];
    cg1 = sc[0][5];
    sa2 = sc[1][0];
    ca2 = sc[1][1];
    sb2 = sc[1][2];
    cb2 = sc[1][3];
    sg2 = sc[1][4];
    cg2 = sc[1][5];
    sa3 = sc[2][0];
    ca3 = sc[2][1];
    sb3 = sc[2][2];
    cb3 = sc[2][3];
    sg3 = sc[2][4];
    cg3 = sc[2][5];
    sa4 = sc[3][0];
    ca4 = sc[3][1];
    sb4 = sc[3][2];
    cb4 = sc[3][3];
    sg4 = sc[3][4];
    cg4 = sc[3][5];

    // compute sines and cosines of theta

    st1 = sct[0];
    ct1 = sct[1];
    st2 = sct[2];
    ct2 = sct[3];
    st3 = sct[4];
    ct3 = sct[5];
    st4 = sct[6];
    ct4 = sct[7];

    // st1 = sinf(theta[0]);
    // ct1 = cosf(theta[0]);
    // st2 = sinf(theta[1]);
    // ct2 = cosf(theta[1]);
    // st3 = sinf(theta[2]);
    // ct3 = cosf(theta[2]);
    // st4 = sinf(theta[3]);
    // ct4 = cosf(theta[3]);

    //extract mass parameters
    mp = parammass[0];
    // mobile platform % Ip = mass(2);
    // not needed % mj = mass(3);
    // of the joint not needed
    mb = parammass[3];
    Ib = parammass[4];
    ml = parammass[5];

    // extract this
    //float l = (*params)[i][4];
    float b = (*params)[0][4];
    for (int i = 0; i < 7; i++)
    {
        (*coriolis_vector)[i] = 0.0;
        (*gravity_vector)[i] = 0.0;
        for (int j = 0; j < 7; j++)
            (*mass_matrix)[i][j] = 0.0;
    }
    // compute the mass matrix

    (*mass_matrix)[0][0] = Ib + 
                           (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) / 4;
    // (*mass_matrix)[0][1] = 0;
    // (*mass_matrix)[0][2] = 0;
    // (*mass_matrix)[0][3] = 0;
    (*mass_matrix)[0][4] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4;
    (*mass_matrix)[0][5] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4;
    (*mass_matrix)[0][6] = (b * ml * sb1 * st1) / 4;

    // (*mass_matrix)[1][0] = 0;
    (*mass_matrix)[1][1] = Ib + 
                           (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) / 4;
    // (*mass_matrix)[1][2] = 0;
    // (*mass_matrix)[1][3] = 0;
    (*mass_matrix)[1][4] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4;
    (*mass_matrix)[1][5] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4;
    (*mass_matrix)[1][6] = (b * ml * sb2 * st2) / 4;

    // (*mass_matrix)[2][0] = 0;
    // (*mass_matrix)[2][1] = 0;
    (*mass_matrix)[2][2] = Ib + 
                           (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) / 4;
    // (*mass_matrix)[2][3] = 0;
    (*mass_matrix)[2][4] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4;
    (*mass_matrix)[2][5] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4;
    (*mass_matrix)[2][6] = (b * ml * sb3 * st3) / 4;

    // (*mass_matrix)[3][0] = 0;
    // (*mass_matrix)[3][1] = 0;
    // (*mass_matrix)[3][2] = 0;
    (*mass_matrix)[3][3] = Ib + 
                           (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) / 4;
    (*mass_matrix)[3][4] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4;
    (*mass_matrix)[3][5] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4;
    (*mass_matrix)[3][6] = (b * ml * sb4 * st4) / 4;

    (*mass_matrix)[4][0] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4;
    (*mass_matrix)[4][1] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4;
    (*mass_matrix)[4][2] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4;
    (*mass_matrix)[4][3] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4;
    (*mass_matrix)[4][4] = ml + mp;
    // (*mass_matrix)[4][5] = 0;
    // (*mass_matrix)[4][6] = 0;

    (*mass_matrix)[5][0] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4;
    (*mass_matrix)[5][1] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4;
    (*mass_matrix)[5][2] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4;
    (*mass_matrix)[5][3] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4;
    // (*mass_matrix)[5][4] = 0;
    (*mass_matrix)[5][5] = ml + mp;
    // (*mass_matrix)[5][6] = 0;

    (*mass_matrix)[6][0] = (b * ml * sb1 * st1) / 4;
    (*mass_matrix)[6][1] = (b * ml * sb2 * st2) / 4;
    (*mass_matrix)[6][2] = (b * ml * sb3 * st3) / 4;
    (*mass_matrix)[6][3] = (b * ml * sb4 * st4) / 4;
    // (*mass_matrix)[6][4] = 0;
    // (*mass_matrix)[6][5] = 0;
    (*mass_matrix)[6][6] = ml + mp;

    // compute the centrifugal and coriolis vector

    // (*coriolis_vector)[0] = 0.0;
    // (*coriolis_vector)[1] = 0.0;
    // (*coriolis_vector)[2] = 0.0;
    // (*coriolis_vector)[3] = 0.0;
    (*coriolis_vector)[4] = 
        (b * ml * (powf(dtheta[0], 2.0) * sa1 * st1 + 
        powf(dtheta[1], 2.0) * sa2 * st2 + powf(dtheta[2], 2.0) * sa3 * st3 + 
        powf(dtheta[3], 2.0) * sa4 * st4 - powf(dtheta[0], 2.0) * ca1 * cb1 * 
        ct1 - powf(dtheta[1], 2.0) * ca2 * cb2 * ct2 - powf(dtheta[2], 2.0) * 
        ca3 * cb3 * ct3 - powf(dtheta[3], 2.0) * ca4 * cb4 * ct4)) / 4;
    (*coriolis_vector)[5] = -(b * ml * (powf(dtheta[0], 2.0) * ca1 * st1 + 
        powf(dtheta[1], 2.0) * ca2 * st2 + powf(dtheta[2], 2.0) * ca3 * st3 + 
        powf(dtheta[3], 2.0) * ca4 * st4 + powf(dtheta[0], 2.0) * cb1 * sa1 * 
        ct1 + powf(dtheta[1], 2.0) * cb2 * sa2 * ct2 + powf(dtheta[2], 2.0) * 
        cb3 * sa3 * ct3 + powf(dtheta[3], 2.0) * cb4 * sa4 * ct4)) / 4;
    (*coriolis_vector)[6] = (b * ml * (sb1 * ct1 * powf(dtheta[0], 2.0) + 
        sb2 * ct2 * powf(dtheta[1], 2.0) + sb3 * ct3 * powf(dtheta[2], 2.0) + 
        sb4 * ct4 * powf(dtheta[3], 2.0))) / 4;

    (*gravity_vector)[0] = (b * g * mb * sb1 * st1) / 2 + (b * g * ml * sb1 * 
        st1) / 2.0;
    (*gravity_vector)[1] = (b * g * mb * sb2 * st2) / 2 + (b * g * ml * sb2 * 
        st2) / 2.0;
    (*gravity_vector)[2] = (b * g * mb * sb3 * st3) / 2 + (b * g * ml * sb3 *
        st3) / 2.0;
    (*gravity_vector)[3] = (b * g * mb * sb4 * st4) / 2 + (b * g * ml * sb4 * 
        st4) / 2.0;
    // (*gravity_vector)[4] = 0;
    // (*gravity_vector)[5] = 0;
    (*gravity_vector)[6] = g * (2 * ml + mp);
}

void ragnarMassCentrGravf(
    float theta[4], float dtheta[4], float (*params)[4][8], float parammass[6],
    float sc[4][6], float (*mass_matrix)[7][7], float (*coriolis_vector)[7], 
    float (*gravity_vector)[7], float sct[8], float g)
{
    float sa1, ca1, sb1, cb1, sg1, cg1, sa2, ca2, sb2, cb2, sg2, cg2;
    float sa3, ca3, sb3, cb3, sg3, cg3, sa4, ca4, sb4, cb4, sg4, cg4;
    float st1, ct1, st2, ct2, st3, ct3, st4, ct4;
    float mp, mb, Ib, ml;
    // extract precomputed sines and cosines
    sa1 = sc[0][0];
    ca1 = sc[0][1];
    sb1 = sc[0][2];
    cb1 = sc[0][3];
    sg1 = sc[0][4];
    cg1 = sc[0][5];
    sa2 = sc[1][0];
    ca2 = sc[1][1];
    sb2 = sc[1][2];
    cb2 = sc[1][3];
    sg2 = sc[1][4];
    cg2 = sc[1][5];
    sa3 = sc[2][0];
    ca3 = sc[2][1];
    sb3 = sc[2][2];
    cb3 = sc[2][3];
    sg3 = sc[2][4];
    cg3 = sc[2][5];
    sa4 = sc[3][0];
    ca4 = sc[3][1];
    sb4 = sc[3][2];
    cb4 = sc[3][3];
    sg4 = sc[3][4];
    cg4 = sc[3][5];

    // compute sines and cosines of theta

    st1 = sct[0];
    ct1 = sct[1];
    st2 = sct[2];
    ct2 = sct[3];
    st3 = sct[4];
    ct3 = sct[5];
    st4 = sct[6];
    ct4 = sct[7];

    // st1 = sinf(theta[0]);
    // ct1 = cosf(theta[0]);
    // st2 = sinf(theta[1]);
    // ct2 = cosf(theta[1]);
    // st3 = sinf(theta[2]);
    // ct3 = cosf(theta[2]);
    // st4 = sinf(theta[3]);
    // ct4 = cosf(theta[3]);

    //extract mass parameters
    mp = parammass[0];
    // mobile platform % Ip = mass(2);
    // not needed % mj = mass(3);
    // of the joint not needed
    mb = parammass[3];
    Ib = parammass[4];
    ml = parammass[5];

    // extract this
    //float l = (*params)[i][4];
    float b = (*params)[0][4];
    for (int i = 0; i < 7; i++)
    {
        (*coriolis_vector)[i] = 0.0;
        (*gravity_vector)[i] = 0.0;
        for (int j = 0; j < 7; j++)
            (*mass_matrix)[i][j] = 0.0;
    }
    // compute the mass matrix

    (*mass_matrix)[0][0] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // (*mass_matrix)[0][1] = 0;
    // (*mass_matrix)[0][2] = 0;
    // (*mass_matrix)[0][3] = 0;
    (*mass_matrix)[0][4] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4.0;
    (*mass_matrix)[0][5] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4.0;
    (*mass_matrix)[0][6] = (b * ml * sb1 * st1) / 4;

    // (*mass_matrix)[1][0] = 0;
    (*mass_matrix)[1][1] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // (*mass_matrix)[1][2] = 0;
    // (*mass_matrix)[1][3] = 0;
    (*mass_matrix)[1][4] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4;
    (*mass_matrix)[1][5] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4;
    (*mass_matrix)[1][6] = (b * ml * sb2 * st2) / 4;

    // (*mass_matrix)[2][0] = 0;
    // (*mass_matrix)[2][1] = 0;
    (*mass_matrix)[2][2] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // (*mass_matrix)[2][3] = 0;
    (*mass_matrix)[2][4] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4;
    (*mass_matrix)[2][5] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4;
    (*mass_matrix)[2][6] = (b * ml * sb3 * st3) / 4;

    // (*mass_matrix)[3][0] = 0;
    // (*mass_matrix)[3][1] = 0;
    // (*mass_matrix)[3][2] = 0;
    (*mass_matrix)[3][3] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    (*mass_matrix)[3][4] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4;
    (*mass_matrix)[3][5] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4;
    (*mass_matrix)[3][6] = (b * ml * sb4 * st4) / 4;

    (*mass_matrix)[4][0] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4;
    (*mass_matrix)[4][1] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4;
    (*mass_matrix)[4][2] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4;
    (*mass_matrix)[4][3] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4;
    (*mass_matrix)[4][4] = ml + mp;
    // (*mass_matrix)[4][5] = 0;
    // (*mass_matrix)[4][6] = 0;

    (*mass_matrix)[5][0] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4;
    (*mass_matrix)[5][1] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4;
    (*mass_matrix)[5][2] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4;
    (*mass_matrix)[5][3] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4;
    // (*mass_matrix)[5][4] = 0;
    (*mass_matrix)[5][5] = ml + mp;
    // (*mass_matrix)[5][6] = 0;

    (*mass_matrix)[6][0] = (b * ml * sb1 * st1) / 4;
    (*mass_matrix)[6][1] = (b * ml * sb2 * st2) / 4;
    (*mass_matrix)[6][2] = (b * ml * sb3 * st3) / 4;
    (*mass_matrix)[6][3] = (b * ml * sb4 * st4) / 4;
    // (*mass_matrix)[6][4] = 0;
    // (*mass_matrix)[6][5] = 0;
    (*mass_matrix)[6][6] = ml + mp;

    // compute the centrifugal and coriolis vector

    // (*coriolis_vector)[0] = 0.0;
    // (*coriolis_vector)[1] = 0.0;
    // (*coriolis_vector)[2] = 0.0;
    // (*coriolis_vector)[3] = 0.0;
    (*coriolis_vector)[4] = (b * ml * (powf(dtheta[0], 2.0) * sa1 * st1 + 
        powf(dtheta[1], 2.0) * sa2 * st2 + powf(dtheta[2], 2.0) * sa3 * st3 + 
        powf(dtheta[3], 2.0) * sa4 * st4 - powf(dtheta[0], 2.0) * ca1 * cb1 * 
        ct1 - powf(dtheta[1], 2.0) * ca2 * cb2 * ct2 - powf(dtheta[2], 2.0) * 
        ca3 * cb3 * ct3 - powf(dtheta[3], 2.0) * ca4 * cb4 * ct4)) / 4.0;
    (*coriolis_vector)[5] = -(b * ml * (powf(dtheta[0], 2.0) * ca1 * st1 + 
        powf(dtheta[1], 2.0) * ca2 * st2 + powf(dtheta[2], 2.0) * ca3 * st3 + 
        powf(dtheta[3], 2.0) * ca4 * st4 + powf(dtheta[0], 2.0) * cb1 * sa1 * 
        ct1 + powf(dtheta[1], 2.0) * cb2 * sa2 * ct2 + powf(dtheta[2], 2.0) * 
        cb3 * sa3 * ct3 + powf(dtheta[3], 2.0) * cb4 * sa4 * ct4)) / 4.0;
    (*coriolis_vector)[6] = (b * ml * (sb1 * ct1 * powf(dtheta[0], 2.0) + sb2 * 
        ct2 * powf(dtheta[1], 2.0) + sb3 * ct3 * powf(dtheta[2], 2.0) + sb4 * 
        ct4 * powf(dtheta[3], 2.0))) / 4.0;

    (*gravity_vector)[0] = (b * g * mb * sb1 * st1) / 2 + (b * g * ml * sb1 * 
        st1) / 2.0;
    (*gravity_vector)[1] = (b * g * mb * sb2 * st2) / 2 + (b * g * ml * sb2 * 
        st2) / 2.0;
    (*gravity_vector)[2] = (b * g * mb * sb3 * st3) / 2 + (b * g * ml * sb3 * 
        st3) / 2.0;
    (*gravity_vector)[3] = (b * g * mb * sb4 * st4) / 2 + (b * g * ml * sb4 * 
        st4) / 2.0;
    // (*gravity_vector)[4] = 0;
    // (*gravity_vector)[5] = 0;
    (*gravity_vector)[6] = g * (2 * ml + mp);
}

void ragnarMassvCentrGrav(
    double theta[4], double dtheta[4], double ddq[7], float params[4][8], 
    double parammass[6], double sc[4][6], double (*mass_vector)[7], 
    double (*coriolis_vector)[7], double (*gravity_vector)[7], double sct[8],
    double g)
{
    double sa1, ca1, sb1, cb1, sg1, cg1, sa2, ca2, sb2, cb2, sg2, cg2;
    double sa3, ca3, sb3, cb3, sg3, cg3, sa4, ca4, sb4, cb4, sg4, cg4;
    double st1, ct1, st2, ct2, st3, ct3, st4, ct4;
    double mp, mb, Ib, ml;
    // extract precomputed sines and cosines
    sa1 = sc[0][0];
    ca1 = sc[0][1];
    sb1 = sc[0][2];
    cb1 = sc[0][3];
    sg1 = sc[0][4];
    cg1 = sc[0][5];
    sa2 = sc[1][0];
    ca2 = sc[1][1];
    sb2 = sc[1][2];
    cb2 = sc[1][3];
    sg2 = sc[1][4];
    cg2 = sc[1][5];
    sa3 = sc[2][0];
    ca3 = sc[2][1];
    sb3 = sc[2][2];
    cb3 = sc[2][3];
    sg3 = sc[2][4];
    cg3 = sc[2][5];
    sa4 = sc[3][0];
    ca4 = sc[3][1];
    sb4 = sc[3][2];
    cb4 = sc[3][3];
    sg4 = sc[3][4];
    cg4 = sc[3][5];

    // compute sines and cosines of theta

    st1 = sct[0];
    ct1 = sct[1];
    st2 = sct[2];
    ct2 = sct[3];
    st3 = sct[4];
    ct3 = sct[5];
    st4 = sct[6];
    ct4 = sct[7];

    // st1 = sinf(theta[0]);
    // ct1 = cosf(theta[0]);
    // st2 = sinf(theta[1]);
    // ct2 = cosf(theta[1]);
    // st3 = sinf(theta[2]);
    // ct3 = cosf(theta[2]);
    // st4 = sinf(theta[3]);
    // ct4 = cosf(theta[3]);

    //extract mass parameters
    mp = parammass[0];
    // mobile platform % Ip = mass(2);
    // not needed % mj = mass(3);
    // of the joint not needed
    mb = parammass[3];
    Ib = parammass[4];
    ml = parammass[5];
    double mass_matrix[7][7];
    // extract this
    //float l = (*params)[i][4];
    float b = params[0][4];
    for (int i = 0; i < 7; i++)
    {
        (*coriolis_vector)[i] = 0.0;
        (*gravity_vector)[i] = 0.0;
        for (int j = 0; j < 7; j++)
            mass_matrix[i][j] = 0.0;
    }
    // compute the mass matrix

    mass_matrix[0][0] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[0][1] = 0;
    // mass_matrix[0][2] = 0;
    // mass_matrix[0][3] = 0;
    mass_matrix[0][4] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4.0;
    mass_matrix[0][5] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4.0;
    mass_matrix[0][6] = (b * ml * sb1 * st1) / 4.0;

    // mass_matrix[1][0] = 0;
    mass_matrix[1][1] = Ib + (powf(b, 2.0) * mb) / 4.0 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[1][2] = 0;
    // mass_matrix[1][3] = 0;
    mass_matrix[1][4] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4.0;
    mass_matrix[1][5] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4.0;
    mass_matrix[1][6] = (b * ml * sb2 * st2) / 4.0;

    // mass_matrix[2][0] = 0;
    // mass_matrix[2][1] = 0;
    mass_matrix[2][2] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[2][3] = 0;
    mass_matrix[2][4] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4.0;
    mass_matrix[2][5] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4.0;
    mass_matrix[2][6] = (b * ml * sb3 * st3) / 4.0;

    // mass_matrix[3][0] = 0;
    // mass_matrix[3][1] = 0;
    // mass_matrix[3][2] = 0;
    mass_matrix[3][3] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    mass_matrix[3][4] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4.0;
    mass_matrix[3][5] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4.0;
    mass_matrix[3][6] = (b * ml * sb4 * st4) / 4;

    mass_matrix[4][0] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4.0;
    mass_matrix[4][1] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4.0;
    mass_matrix[4][2] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4.0;
    mass_matrix[4][3] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4.0;
    mass_matrix[4][4] = ml + mp;
    // mass_matrix[4][5] = 0;
    // mass_matrix[4][6] = 0;

    mass_matrix[5][0] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4.0;
    mass_matrix[5][1] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4.0;
    mass_matrix[5][2] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4.0;
    mass_matrix[5][3] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4.0;
    // mass_matrix[5][4] = 0;
    mass_matrix[5][5] = ml + mp;
    // mass_matrix[5][6] = 0;

    mass_matrix[6][0] = (b * ml * sb1 * st1) / 4.0;
    mass_matrix[6][1] = (b * ml * sb2 * st2) / 4.0;
    mass_matrix[6][2] = (b * ml * sb3 * st3) / 4.0;
    mass_matrix[6][3] = (b * ml * sb4 * st4) / 4.0;
    // mass_matrix[6][4] = 0;
    // mass_matrix[6][5] = 0;
    mass_matrix[6][6] = ml + mp;

    // Compute M * ddq
    double temp = 0.0;
    for (int i = 0; i < 7; i++)
    {
        temp = 0.0;
        for (int j = 0; j < 7; j++)
        {
            temp = temp + mass_matrix[i][j] * ddq[j];
        }
        (*mass_vector)[i] = temp;
    }

    // compute the centrifugal and coriolis vector

    // (*coriolis_vector)[0] = 0.0;
    // (*coriolis_vector)[1] = 0.0;
    // (*coriolis_vector)[2] = 0.0;
    // (*coriolis_vector)[3] = 0.0;
    (*coriolis_vector)[4] = (b * ml * (powf(dtheta[0], 2.0) * sa1 * st1 + 
        powf(dtheta[1], 2.0) * sa2 * st2 + powf(dtheta[2], 2.0) * sa3 * st3 + 
        powf(dtheta[3], 2.0) * sa4 * st4 - powf(dtheta[0], 2.0) * ca1 * cb1 * 
        ct1 - powf(dtheta[1], 2.0) * ca2 * cb2 * ct2 - powf(dtheta[2], 2.0) * 
        ca3 * cb3 * ct3 - powf(dtheta[3], 2.0) * ca4 * cb4 * ct4)) / 4.0;
    (*coriolis_vector)[5] = -(b * ml * (powf(dtheta[0], 2.0) * ca1 * st1 + 
        powf(dtheta[1], 2.0) * ca2 * st2 + powf(dtheta[2], 2.0) * ca3 * st3 + 
        powf(dtheta[3], 2.0) * ca4 * st4 + powf(dtheta[0], 2.0) * cb1 * sa1 * 
        ct1 + powf(dtheta[1], 2.0) * cb2 * sa2 * ct2 + powf(dtheta[2], 2.0) * 
        cb3 * sa3 * ct3 + powf(dtheta[3], 2.0) * cb4 * sa4 * ct4)) / 4.0;
    (*coriolis_vector)[6] = (b * ml * (sb1 * ct1 * powf(dtheta[0], 2.0) + sb2 * 
        ct2 * powf(dtheta[1], 2.0) + sb3 * ct3 * powf(dtheta[2], 2.0) + sb4 * 
        ct4 * powf(dtheta[3], 2.0))) / 4.0;

    (*gravity_vector)[0] = (b * g * mb * sb1 * st1) / 2.0 + (b * g * ml * sb1 *
        st1) / 2.0;
    (*gravity_vector)[1] = (b * g * mb * sb2 * st2) / 2.0 + (b * g * ml * sb2 * 
        st2) / 2.0;
    (*gravity_vector)[2] = (b * g * mb * sb3 * st3) / 2.0 + (b * g * ml * sb3 * 
        st3) / 2.0;
    (*gravity_vector)[3] = (b * g * mb * sb4 * st4) / 2.0 + (b * g * ml * sb4 * 
        st4) / 2.0;
    // (*gravity_vector)[4] = 0;
    // (*gravity_vector)[5] = 0;
    (*gravity_vector)[6] = g * (2 * ml + mp);
}

void ragnarMassvCentrGravf(
    float theta[4], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6], float (*mass_vector)[7], 
    float (*coriolis_vector)[7], float (*gravity_vector)[7], float sct[8],
    float g)
{
    float sa1, ca1, sb1, cb1, sg1, cg1, sa2, ca2, sb2, cb2, sg2, cg2;
    float sa3, ca3, sb3, cb3, sg3, cg3, sa4, ca4, sb4, cb4, sg4, cg4;
    float st1, ct1, st2, ct2, st3, ct3, st4, ct4;
    float mp, mb, Ib, ml;
    // extract precomputed sines and cosines
    sa1 = sc[0][0];
    ca1 = sc[0][1];
    sb1 = sc[0][2];
    cb1 = sc[0][3];
    sg1 = sc[0][4];
    cg1 = sc[0][5];
    sa2 = sc[1][0];
    ca2 = sc[1][1];
    sb2 = sc[1][2];
    cb2 = sc[1][3];
    sg2 = sc[1][4];
    cg2 = sc[1][5];
    sa3 = sc[2][0];
    ca3 = sc[2][1];
    sb3 = sc[2][2];
    cb3 = sc[2][3];
    sg3 = sc[2][4];
    cg3 = sc[2][5];
    sa4 = sc[3][0];
    ca4 = sc[3][1];
    sb4 = sc[3][2];
    cb4 = sc[3][3];
    sg4 = sc[3][4];
    cg4 = sc[3][5];

    // compute sines and cosines of theta

    st1 = sct[0];
    ct1 = sct[1];
    st2 = sct[2];
    ct2 = sct[3];
    st3 = sct[4];
    ct3 = sct[5];
    st4 = sct[6];
    ct4 = sct[7];

    // st1 = sinf(theta[0]);
    // ct1 = cosf(theta[0]);
    // st2 = sinf(theta[1]);
    // ct2 = cosf(theta[1]);
    // st3 = sinf(theta[2]);
    // ct3 = cosf(theta[2]);
    // st4 = sinf(theta[3]);
    // ct4 = cosf(theta[3]);

    //extract mass parameters
    mp = parammass[0];
    // mobile platform % Ip = mass(2);
    // not needed % mj = mass(3);
    // of the joint not needed
    mb = parammass[3];
    Ib = parammass[4];
    ml = parammass[5];
    float mass_matrix[7][7];
    // extract this
    //float l = (*params)[i][4];
    float b = params[0][4];
    for (int i = 0; i < 7; i++)
    {
        (*coriolis_vector)[i] = 0.0;
        (*gravity_vector)[i] = 0.0;
        for (int j = 0; j < 7; j++)
            mass_matrix[i][j] = 0.0;
    }
    // compute the mass matrix

    mass_matrix[0][0] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[0][1] = 0;
    // mass_matrix[0][2] = 0;
    // mass_matrix[0][3] = 0;
    mass_matrix[0][4] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) 
        / 4.0;
    mass_matrix[0][5] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) 
        / 4.0;
    mass_matrix[0][6] = (b * ml * sb1 * st1) / 4.0;

    // mass_matrix[1][0] = 0;
    mass_matrix[1][1] = Ib + (powf(b, 2.0) * mb) / 4.0 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[1][2] = 0;
    // mass_matrix[1][3] = 0;
    mass_matrix[1][4] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4.0;
    mass_matrix[1][5] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4.0;
    mass_matrix[1][6] = (b * ml * sb2 * st2) / 4.0;

    // mass_matrix[2][0] = 0;
    // mass_matrix[2][1] = 0;
    mass_matrix[2][2] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    // mass_matrix[2][3] = 0;
    mass_matrix[2][4] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4.0;
    mass_matrix[2][5] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4.0;
    mass_matrix[2][6] = (b * ml * sb3 * st3) / 4.0;

    // mass_matrix[3][0] = 0;
    // mass_matrix[3][1] = 0;
    // mass_matrix[3][2] = 0;
    mass_matrix[3][3] = Ib + (powf(b, 2.0) * mb) / 4 + (powf(b, 2.0) * ml) 
        / 4.0;
    mass_matrix[3][4] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4.0;
    mass_matrix[3][5] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4.0;
    mass_matrix[3][6] = (b * ml * sb4 * st4) / 4;

    mass_matrix[4][0] = -(b * ml * (sa1 * ct1 + ca1 * cb1 * st1)) / 4.0;
    mass_matrix[4][1] = -(b * ml * (sa2 * ct2 + ca2 * cb2 * st2)) / 4.0;
    mass_matrix[4][2] = -(b * ml * (sa3 * ct3 + ca3 * cb3 * st3)) / 4.0;
    mass_matrix[4][3] = -(b * ml * (sa4 * ct4 + ca4 * cb4 * st4)) / 4.0;
    mass_matrix[4][4] = ml + mp;
    // mass_matrix[4][5] = 0;
    // mass_matrix[4][6] = 0;

    mass_matrix[5][0] = (b * ml * (ca1 * ct1 - cb1 * sa1 * st1)) / 4.0;
    mass_matrix[5][1] = (b * ml * (ca2 * ct2 - cb2 * sa2 * st2)) / 4.0;
    mass_matrix[5][2] = (b * ml * (ca3 * ct3 - cb3 * sa3 * st3)) / 4.0;
    mass_matrix[5][3] = (b * ml * (ca4 * ct4 - cb4 * sa4 * st4)) / 4.0;
    // mass_matrix[5][4] = 0;
    mass_matrix[5][5] = ml + mp;
    // mass_matrix[5][6] = 0;

    mass_matrix[6][0] = (b * ml * sb1 * st1) / 4.0;
    mass_matrix[6][1] = (b * ml * sb2 * st2) / 4.0;
    mass_matrix[6][2] = (b * ml * sb3 * st3) / 4.0;
    mass_matrix[6][3] = (b * ml * sb4 * st4) / 4.0;
    // mass_matrix[6][4] = 0;
    // mass_matrix[6][5] = 0;
    mass_matrix[6][6] = ml + mp;

    // Compute M * ddq
    float temp = 0.0;
    for (int i = 0; i < 7; i++)
    {
        temp = 0.0;
        for (int j = 0; j < 7; j++)
        {
            temp = temp + mass_matrix[i][j] * ddq[j];
        }
        (*mass_vector)[i] = temp;
    }

    // compute the centrifugal and coriolis vector

    // (*coriolis_vector)[0] = 0.0;
    // (*coriolis_vector)[1] = 0.0;
    // (*coriolis_vector)[2] = 0.0;
    // (*coriolis_vector)[3] = 0.0;
    (*coriolis_vector)[4] = (b * ml * (powf(dtheta[0], 2.0) * sa1 * st1 + 
        powf(dtheta[1], 2.0) * sa2 * st2 + powf(dtheta[2], 2.0) * sa3 * st3 + 
        powf(dtheta[3], 2.0) * sa4 * st4 - powf(dtheta[0], 2.0) * ca1 * cb1 * 
        ct1 - powf(dtheta[1], 2.0) * ca2 * cb2 * ct2 - powf(dtheta[2], 2.0) * 
        ca3 * cb3 * ct3 - powf(dtheta[3], 2.0) * ca4 * cb4 * ct4)) / 4.0;
    (*coriolis_vector)[5] = -(b * ml * (powf(dtheta[0], 2.0) * ca1 * st1 + 
        powf(dtheta[1], 2.0) * ca2 * st2 + powf(dtheta[2], 2.0) * ca3 * st3 + 
        powf(dtheta[3], 2.0) * ca4 * st4 + powf(dtheta[0], 2.0) * cb1 * sa1 * 
        ct1 + powf(dtheta[1], 2.0) * cb2 * sa2 * ct2 + powf(dtheta[2], 2.0) * 
        cb3 * sa3 * ct3 + powf(dtheta[3], 2.0) * cb4 * sa4 * ct4)) / 4.0;
    (*coriolis_vector)[6] = (b * ml * (sb1 * ct1 * powf(dtheta[0], 2.0) + sb2 * 
        ct2 * powf(dtheta[1], 2.0) + sb3 * ct3 * powf(dtheta[2], 2.0) + sb4 * 
        ct4 * powf(dtheta[3], 2.0))) / 4.0;

    (*gravity_vector)[0] = (b * g * mb * sb1 * st1) / 2.0 + (b * g * ml * sb1 * 
        st1) / 2.0;
    (*gravity_vector)[1] = (b * g * mb * sb2 * st2) / 2.0 + (b * g * ml * sb2 * 
        st2) / 2.0;
    (*gravity_vector)[2] = (b * g * mb * sb3 * st3) / 2.0 + (b * g * ml * sb3 * 
        st3) / 2.0;
    (*gravity_vector)[3] = (b * g * mb * sb4 * st4) / 2.0 + (b * g * ml * sb4 * 
        st4) / 2.0;
    // (*gravity_vector)[4] = 0;
    // (*gravity_vector)[5] = 0;
    (*gravity_vector)[6] = g * (2 * ml + mp);
}

void ragnarTorques(
    double q[7], double dtheta[4], double ddq[7], float params[4][8], 
    double parammass[6], double sc[4][6],double sct[8], double g, 
    double (*torque)[4])
{
    double massvector[7];
    double coriolis[7];
    double grav[7];
    double theta[4];
    theta[0] = q[0];
    theta[1] = q[1];
    theta[2] = q[2];
    theta[3] = q[3];
    // get the mass vector Mq*ddq , coriolis and centrifugal vector and gravity
    // vector
    ragnarMassvCentrGrav(
        theta, dtheta, ddq, params, parammass, sc, &massvector, &coriolis, 
            &grav, sct, g);
    //**** Create the right hand side of the equation
    double f[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 7; i++)
        f[i] = massvector[i] + coriolis[i] + grav[i];

    // get the constraint jacobian
    double A[4][3];
    double B[4][4];

    ragnarAB(q, params, sct, sc, &A, &B);

    //**** Create the left hand side matrix to invert
    double MC[7][8];
    //  Serial.println("CM matrix to show if correct");
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MC[i][j] = 0.0;
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if (j < 4)
            {
                if (i < 4)
                { // create the diagonal matrix inside the matrix
                    if (i == j)
                        MC[i][j] = 1;
                }
            }
            else
            { // put the Cq' here = [B -A]'
                if (i < 4)
                    MC[i][j] = B[j - 4][i]; // transpose exchange rows for columns
                else
                    MC[i][j] = -A[j - 4][i - 4];
            }
        }
    }
    double MCt[8][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MCt[j][i] = MC[i][j];

    //  Serial.println("finish CM" );
    Matrix<double> mm(7, 8, (double *)MC);
    Matrix<double> mmt(8, 7, (double *)MCt);
    // Matrix<double> CMi = Matrix<double>::transpose(mm) * 
    //     Matrix<double>::inv(mm * Matrix<double>::transpose(mm) );
    Matrix<double> temp = mm * mmt;

    // Test to do the inverse with float precision
    float temp_1[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_1[i][j] = float(temp._entity[i][j]);
    
    Matrix<float> temp_2(7, 7, (float *)temp_1);
    Matrix<float> iitemp = Matrix<float>::inv(temp_2);
    // Return to a double Matrix to do the multiplication to the right hand
    // side. 
    double temp_3[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_3[i][j] = double(iitemp._entity[i][j]);

    /// Matrix<double> itemp = Matrix<double>::inv(temp);
    Matrix<double> itemp(7, 7, (double *)temp_3);

    //  Matrix<double> CMi = mmt * Matrix<double>::inv(mm * mmt );
    Matrix<double> CMi = mmt * itemp;
    Matrix<double> fvector(7, 1, (double *)f);
    Matrix<double> tau_lambda = CMi * fvector;

    (*torque)[0] = tau_lambda._entity[0][0];
    (*torque)[1] = tau_lambda._entity[1][0];
    (*torque)[2] = tau_lambda._entity[2][0];
    (*torque)[3] = tau_lambda._entity[3][0];
}
void ragnarTorquesfimp(
    float q[7], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6],float sct[8], float g, 
    float (*torque)[4], float impednce[3])
{
    float massvector[7];
    float coriolis[7];
    float grav[7];
    float theta[4];
    theta[0] = q[0];
    theta[1] = q[1];
    theta[2] = q[2];
    theta[3] = q[3];
    // get the mass vector Mq*ddq , coriolis and centrifugal vector and gravity
    // vector
    ragnarMassvCentrGravf(
        theta, dtheta, ddq, params, parammass, sc, &massvector, &coriolis, 
        &grav, sct, g);
    // TEST ADDING THE FRICTION PARAMETERS 
    //friction = [0.6841, 0.0098; 0.8224, 0.018; 0.6646, 0.0188; 0.6911, 0.0259];

    // float f1f2[2][4] = {{0.6841, 0.8224, 0.6646, 0.6911},
    //                     {0.0098, 0.018, 0.0188, 0.0259}};
    float f1f2[2][4] = {{0.6841, 0.7224, 0.6646, 0.6911},
                        {0.0098, 0.018, 0.0188, 0.0259}};

    float friction[7]; 
    for (int i=0; i<4; i++){
        if (dtheta[i] > 0.0) 
            friction[i] = f1f2[0][i] * (1-expf(-3.0 * abs(dtheta[i])));
        else if (dtheta[i] < 0.0) 
            friction[i] = -f1f2[0][i] * (1-expf(-3.0 * abs(dtheta[i])));
        else friction[i] = 0.0; 
        
        friction[i] = friction[i] + f1f2[1][i]*dtheta[i];
    }
    for (int i=4; i<7; i++) friction[i] = 0; 
    //**** Create the right hand side of the equation
    float f[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    float impd[7]; 
    for (int i=0; i<4; i++) impd[i] = 0.0;
    for (int i=4; i<7; i++) impd[i] = impednce[i-4]; 
    
    for (int i = 0; i < 7; i++)
        f[i] = massvector[i] + coriolis[i] + grav[i] + impd[i]; // + friction[i]; 
    // Matrix<float> massv(7,1, (float*) f);
    // massv.show();
    
    // get the constraint jacobian
    float A[4][3];
    float B[4][4];

    ragnarABf(q, params, sct, sc, &A, &B);

    //**** Create the left hand side matrix to invert
    float MC[7][8];
    //  Serial.println("CM matrix to show if correct");
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MC[i][j] = 0.0;
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if (j < 4)
            {
                if (i < 4)
                { // create the diagonal matrix inside the matrix
                    if (i == j)
                        MC[i][j] = 1;
                }
            }
            else
            { // put the Cq' here = [B -A]'
                if (i < 4)
                    // transpose exchange rows for columns
                    MC[i][j] = B[j - 4][i]; 
                else
                    MC[i][j] = -A[j - 4][i - 4];
            }
        }
    }
    // Matrix<float> massv(7,8, (float*) MC);
    // massv.show();
    
    float MCt[8][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MCt[j][i] = MC[i][j];

    //  Serial.println("finish CM" );
    Matrix<float> mm(7, 8, (float *)MC);
    Matrix<float> mmt(8, 7, (float *)MCt);
    // Matrix<float> CMi = Matrix<float>::transpose(mm) * Matrix<float>::inv(mm * Matrix<float>::transpose(mm) );
    Matrix<float> temp = mm * mmt;

    // Test to do the inverse with float precision
    float temp_1[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_1[i][j] = float(temp._entity[i][j]);
    
    Matrix<float> temp_2(7, 7, (float *)temp_1);
    Matrix<float> iitemp = Matrix<float>::inv(temp_2);
    // Return to a float Matrix to do the multiplication to the right hand
    // side. 
    float temp_3[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_3[i][j] = float(iitemp._entity[i][j]);

    /// Matrix<float> itemp = Matrix<float>::inv(temp);
    Matrix<float> itemp(7, 7, (float *)temp_3);

    //  Matrix<float> CMi = mmt * Matrix<float>::inv(mm * mmt );
    Matrix<float> CMi = mmt * itemp;
    Matrix<float> fvector(7, 1, (float *)f);
    Matrix<float> tau_lambda = CMi * fvector;

    (*torque)[0] = tau_lambda._entity[0][0];
    (*torque)[1] = tau_lambda._entity[1][0];
    (*torque)[2] = tau_lambda._entity[2][0];
    (*torque)[3] = tau_lambda._entity[3][0];

    for(int i=0; i<4; i++) (*torque)[i] = tau_lambda._entity[i][0];
                                           // + friction[i];
}
void ragnarTorquesf(
    float q[7], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6],float sct[8], float g, 
    float (*torque)[4])
{
    float massvector[7];
    float coriolis[7];
    float grav[7];
    float theta[4];
    theta[0] = q[0];
    theta[1] = q[1];
    theta[2] = q[2];
    theta[3] = q[3];
    // get the mass vector Mq*ddq , coriolis and centrifugal vector and gravity
    // vector
    ragnarMassvCentrGravf(
        theta, dtheta, ddq, params, parammass, sc, &massvector, &coriolis, 
        &grav, sct, g);
    // TEST ADDING THE FRICTION PARAMETERS 
    //friction = [0.6841, 0.0098; 0.8224, 0.018; 0.6646, 0.0188; 0.6911, 0.0259];

    // float f1f2[2][4] = {{0.6841, 0.8224, 0.6646, 0.6911},
    //                     {0.0098, 0.018, 0.0188, 0.0259}};
    float f1f2[2][4] = {{0.6841, 0.7224, 0.6646, 0.6911},
                        {0.0098, 0.018, 0.0188, 0.0259}};

    float friction[7]; 
    for (int i=0; i<4; i++){
        if (dtheta[i] > 0.0) 
            friction[i] = f1f2[0][i] * (1-expf(-3.0 * abs(dtheta[i])));
        else if (dtheta[i] < 0.0) 
            friction[i] = -f1f2[0][i] * (1-expf(-3.0 * abs(dtheta[i])));
        else friction[i] = 0.0; 
        
        friction[i] = friction[i] + f1f2[1][i]*dtheta[i];
    }
    for (int i=4; i<7; i++) friction[i] = 0; 
    //**** Create the right hand side of the equation
    float f[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 7; i++)
        f[i] = massvector[i] + coriolis[i] + grav[i]; // + friction[i]; 
    // Matrix<float> massv(7,1, (float*) f);
    // massv.show();
    
    // get the constraint jacobian
    float A[4][3];
    float B[4][4];

    ragnarABf(q, params, sct, sc, &A, &B);

    //**** Create the left hand side matrix to invert
    float MC[7][8];
    //  Serial.println("CM matrix to show if correct");
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MC[i][j] = 0.0;
    for (int i = 0; i < 7; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if (j < 4)
            {
                if (i < 4)
                { // create the diagonal matrix inside the matrix
                    if (i == j)
                        MC[i][j] = 1;
                }
            }
            else
            { // put the Cq' here = [B -A]'
                if (i < 4)
                    // transpose exchange rows for columns
                    MC[i][j] = B[j - 4][i]; 
                else
                    MC[i][j] = -A[j - 4][i - 4];
            }
        }
    }
    // Matrix<float> massv(7,8, (float*) MC);
    // massv.show();
    
    float MCt[8][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 8; j++)
            MCt[j][i] = MC[i][j];

    //  Serial.println("finish CM" );
    Matrix<float> mm(7, 8, (float *)MC);
    Matrix<float> mmt(8, 7, (float *)MCt);
    // Matrix<float> CMi = Matrix<float>::transpose(mm) * Matrix<float>::inv(mm * Matrix<float>::transpose(mm) );
    Matrix<float> temp = mm * mmt;

    // Test to do the inverse with float precision
    float temp_1[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_1[i][j] = float(temp._entity[i][j]);
    
    Matrix<float> temp_2(7, 7, (float *)temp_1);
    Matrix<float> iitemp = Matrix<float>::inv(temp_2);
    // Return to a float Matrix to do the multiplication to the right hand
    // side. 
    float temp_3[7][7];
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 7; j++)
            temp_3[i][j] = float(iitemp._entity[i][j]);

    /// Matrix<float> itemp = Matrix<float>::inv(temp);
    Matrix<float> itemp(7, 7, (float *)temp_3);

    //  Matrix<float> CMi = mmt * Matrix<float>::inv(mm * mmt );
    Matrix<float> CMi = mmt * itemp;
    Matrix<float> fvector(7, 1, (float *)f);
    Matrix<float> tau_lambda = CMi * fvector;

    (*torque)[0] = tau_lambda._entity[0][0];
    (*torque)[1] = tau_lambda._entity[1][0];
    (*torque)[2] = tau_lambda._entity[2][0];
    (*torque)[3] = tau_lambda._entity[3][0];

    for(int i=0; i<4; i++) (*torque)[i] = tau_lambda._entity[i][0];
                                           // + friction[i];
}

bool ragnarLegMassCentrGrav(
    float dtheta[4], float detazeta[4][2], float params[4][8], 
    float parammass[6], float sc[4][6], float sct[8],float g,
    float scez[4][4], int leg_num,float (*mass_matrix)[3][3],
    float (*coriolis_vector)[3], float (*gravity_vector)[3]) 
    {
    float sa = sc[leg_num][0];
    float ca = sc[leg_num][1];
    float sb = sc[leg_num][2];
    float cb = sc[leg_num][3];
    float b = float(params[leg_num][4]); // change notation here, l proximal link, L distal link 
    float L = float(params[leg_num][5]);
    float st = sct[(2*leg_num)];
    float ct = sct[(2*leg_num)+1];
    float sz = scez[leg_num][0];
    float cz = scez[leg_num][1];
    float se = scez[leg_num][2];
    float ce = scez[leg_num][3];
    float dth = dtheta[leg_num];
    float dzeta = detazeta[leg_num][0];
    float deta = detazeta[leg_num][1];
    
    // mobile platform % Ip = mass(2);
    // not needed % mj = mass(3);
    // of the joint not needed
    float mb, Ib, ml;
    mb = parammass[3];
    Ib = parammass[4];
    ml = parammass[5];

    // extract this
    //float l = (*params)[i][4];
    // float b = params[0][4];

    (*mass_matrix)[0][0] = Ib + powf(b,2.0)*ml + (powf(b,2.0)*mb)/4.0 + 
        (powf(L,2.0)*ml*powf(ce,2.0))/4.0 + L*b*ml*ce*cz; 
    (*mass_matrix)[0][1] = (L*ml*ce*(L*ce + 2*b*cz))/4.0;
    (*mass_matrix)[0][2] = -(L*b*ml*se*sz)/2.0;
    (*mass_matrix)[1][0] = (L*ml*ce*(L*ce + 2*b*cz))/4.0;
    (*mass_matrix)[1][1] = (powf(L,2.0)*ml*powf(ce,2.0))/4.0;
    (*mass_matrix)[1][2] = 0.0;
    (*mass_matrix)[2][0] = -(L*b*ml*se*sz)/2.0;
    (*mass_matrix)[2][1] = 0.0;
    (*mass_matrix)[2][2] = (powf(L,2.0)*ml)/4.0;

    (*gravity_vector)[0] = g*ml*(b*sb*st + (L*ce*(sb*ct*sz + sb*cz*st))/2.0) + 
        (b*g*mb*sb*st)/2.0;
    (*gravity_vector)[1] = (L*g*ml*(st*cz+ct*sz)*ce*sb)/2.0;
    (*gravity_vector)[2] = - (L*g*ml*cb*ce)/2.0 - (L*g*ml*se*(sb*st*sz - 
        sb*ct*cz))/2.0;

    (*coriolis_vector)[0] = -(L*ml*(b*powf(deta,2.0)*ce*sz + b*powf(dzeta,2.0)*
        ce*sz + (L*deta*dth*2*se*ce)/2 + (L*deta*dzeta*2*se*ce)/2.0 + 2*b*deta*
        dth*se*cz + 2*b*deta*dzeta*se*cz + 2*b*dth*dzeta*ce*sz))/2.0;
    (*coriolis_vector)[1] = -(L*ml*ce*(- b*sz*powf(dth,2.0) + L*deta*se*dth + 
        L*deta*dzeta*se))/2.0;
    (*coriolis_vector)[2] = (L*ml*se*(L*powf(dth,2.0)*ce + L*powf(dzeta,2.0)*ce 
        + 2*b*powf(dth,2.0)*cz + 2*L*dth*dzeta*ce))/4.0;
}

void ragnarTorquesNf(
    // This is to get the torques with the task space formulation 
    float q[7], float passivef[4][2], float dtheta[4], float ddx[3], 
    float parameter[4][8], float parammass[6], float sc[4][6],float sct[8], 
    float scez[4][4],float g, float (*torque)[4])
{
    float thetaf[4], iJ1[3][3], iJ2[3][3], iJ3[3][3], iJ4[3][3];
    float iJt1[3][3], iJt2[3][3], iJt3[3][3], iJt4[3][3];
    float pJ1[2][3], pJ2[2][3], pJ3[2][3], pJ4[2][3];
    float Af[4][3], Bf[4][4];
    float detazeta[4][2];
    float m1[3][3], m2[3][3], m3[3][3], m4[3][3], c1[3], c2[3], c3[3], c4[3];
    float g1[3], g2[3], g3[3], g4[3];
    // center piece mass matrix 
    float Mp[3][3] = 
        {{parammass[0], 0, 0}, {0, parammass[0], 0}, {0, 0, parammass[0]}};
    // center piece gravity vector 
    float gp[3] = {0, 0, float(g)*parammass[0]};

    thetaf[0] = q[0];
    thetaf[1] = q[1]; 
    thetaf[2] = q[2];
    thetaf[3] = q[3];
    // Compute the inverse Jacobian 
    ragnariJpassivef(parameter, thetaf, passivef, sc, sct, scez, &iJ1, 0);
    ragnariJpassivef(parameter, thetaf, passivef, sc, sct, scez, &iJ2, 1);
    ragnariJpassivef(parameter, thetaf, passivef, sc, sct, scez, &iJ3, 2);
    ragnariJpassivef(parameter, thetaf, passivef, sc, sct, scez, &iJ4, 3);
    // Compute the tranpose of the inverse Jacobian for each arm 
    Atraf(iJ1, &iJt1);
    Atraf(iJ2, &iJt2);
    Atraf(iJ3, &iJt3);
    Atraf(iJ4, &iJt4);
    // Create matrices 
    Matrix<float> miJt1(3, 3, (float*)iJt1);
    Matrix<float> miJt2(3, 3, (float*)iJt2);
    Matrix<float> miJt3(3, 3, (float*)iJt3);
    Matrix<float> miJt4(3, 3, (float*)iJt4);

    Matrix<float> miJ1(3, 3, (float*)iJ1);
    Matrix<float> miJ2(3, 3, (float*)iJ2);
    Matrix<float> miJ3(3, 3, (float*)iJ3);
    Matrix<float> miJ4(3, 3, (float*)iJ4);

    // This function let compute the velocity of the passive joints through a
    // passive Jacobian Jp .- see explanation on the kinematics library 
 
    // Compute the general Jacobian 
    ragnarABf(q, parameter, sct, sc, &Af, &Bf);
    ragnarpassiveJ(Af, Bf, iJ1, &pJ1, 0);
    ragnarpassiveJ(Af, Bf, iJ2, &pJ2, 0);
    ragnarpassiveJ(Af, Bf, iJ3, &pJ3, 1);
    ragnarpassiveJ(Af, Bf, iJ4, &pJ4, 1);
    // Create matrices  
    Matrix<float> mpJ1(2, 3, (float*)pJ1);
    Matrix<float> mpJ2(2, 3, (float*)pJ2);
    Matrix<float> mpJ3(2, 3, (float*)pJ3);
    Matrix<float> mpJ4(2, 3, (float*)pJ4);

    Matrix<float> tempdtheta(4,1,(float*) dtheta);
    // Separate the velocities 
    float dthetaf_1[3] = {tempdtheta._entity[0][0], tempdtheta._entity[1][0], 
        tempdtheta._entity[2][0]};
    float dthetaf_2[3] = {tempdtheta._entity[1][0], tempdtheta._entity[2][0], 
        tempdtheta._entity[3][0]};
    Matrix<float> tempdth1(3, 1, (float*)dthetaf_1);
    Matrix<float> tempdth2(3, 1, (float*)dthetaf_2);
    // Compute the passive velocities 
    Matrix<float> dq1 = mpJ1 * tempdth1;
    Matrix<float> dq2 = mpJ2 * tempdth1;
    Matrix<float> dq3 = mpJ3 * tempdth2;
    Matrix<float> dq4 = mpJ4 * tempdth2;
    // Fill detazeta 
    detazeta[0][0] = dq1._entity[0][0];
    detazeta[0][1] = dq1._entity[1][0];
    detazeta[1][0] = dq2._entity[0][0];
    detazeta[1][1] = dq2._entity[1][0];
    detazeta[2][0] = dq3._entity[0][0];
    detazeta[2][1] = dq3._entity[1][0];
    detazeta[3][0] = dq4._entity[0][0];
    detazeta[3][1] = dq4._entity[1][0];

    // Compute dx 
    Matrix<float> tempAf(4, 3, (float*)Af);
    Matrix<float> tempBf(4, 4, (float*)Bf);
    // get the pseudo inverse of A
    // J = pinv(A)*B
    // A dx = B d\theta
    // dx = pinv(A) * B d\theta
    Matrix<float> at = Matrix<float>::transpose(tempAf);
    Matrix<float> tempaat = at * tempAf;
    Matrix<float> itempaat = Matrix<float>::inv(tempaat);
    Matrix<float> mdx = itempaat * at * tempBf * tempdtheta;

    // Get the mass matrix, coriolis and gravity vectors 
    ragnarLegMassCentrGrav(dtheta, detazeta, parameter, parammass, sc, sct, g, 
        scez, 0, &m1, &c1, &g1);
    ragnarLegMassCentrGrav(dtheta, detazeta, parameter, parammass, sc, sct, g, 
        scez, 1, &m2, &c2, &g2);
    ragnarLegMassCentrGrav(dtheta, detazeta, parameter, parammass, sc, sct, g, 
        scez, 2, &m3, &c3, &g3);
    ragnarLegMassCentrGrav(dtheta, detazeta, parameter, parammass, sc, sct, g, 
        scez, 3, &m4, &c4, &g4);
    // Create matrices 
    Matrix<float> mm1(3, 3, (float*)m1);
    Matrix<float> mm2(3, 3, (float*)m2);
    Matrix<float> mm3(3, 3, (float*)m3);
    Matrix<float> mm4(3, 3, (float*)m4);
    Matrix<float> mc1(3, 1, (float*)c1);
    Matrix<float> mc2(3, 1, (float*)c2);
    Matrix<float> mc3(3, 1, (float*)c3);
    Matrix<float> mc4(3, 1, (float*)c4);
    Matrix<float> mg1(3, 1, (float*)g1);
    Matrix<float> mg2(3, 1, (float*)g2);
    Matrix<float> mg3(3, 1, (float*)g3);
    Matrix<float> mg4(3, 1, (float*)g4);
    // center piece mass matrix 
    Matrix<float> mMp(3, 3, (float*)Mp);
    // Compute the full mass matrix 
    // terms inv(J_i')*m_i*inv(J_i) corresponds to the limb mass matrix in task
    // space 
    Matrix<float> mMM = mMp + miJt1 * mm1 * miJ1 + miJt2 * mm2 * miJ2 
        + miJt3 * mm3 * miJ3 + miJt4 * mm4 * miJ4;
    // center piece gravity vector 
    Matrix<float> mgp(3, 1, (float*)gp);
    // Arms dynamics in task space 
    Matrix<float> arms_dynamics = miJt1 * (mc1 + mg1 ) + miJt2 * (mc2 + mg2)
                                + miJt3 * (mc3 + mg3 ) + miJt4 * (mc4 + mg4);
    // This is a inverse jacobian transpose for the torque computation 
    float IJJT[3][4] = {
        {iJt1[0][0], iJt2[0][0], iJt3[0][0], iJt4[0][0]},
        {iJt1[1][0], iJt2[1][0], iJt3[1][0], iJt4[1][0]},
        {iJt1[2][0], iJt2[2][0], iJt3[2][0], iJt4[2][0]}
        };
    // Create matrix from it 
    Matrix<float> mijjt(3, 4, (float*)IJJT);
    // Create the ddx vector 
    Matrix<float> mddx(3, 1, (float*) ddx);
    //compute right hand of the equation 
    Matrix<float> right_hand = mMM * mddx + mgp + arms_dynamics;

    // Compute the pseudoinverse 
    Matrix<float> mijjtt = Matrix<float>::transpose(mijjt);
    Matrix<float> iterm = Matrix<float>::inv(mijjt * mijjtt);
    Matrix<float> pseudoinverse = mijjtt * iterm;

    // get the torques 
    Matrix<float> torques = pseudoinverse * right_hand;

    // Fill the torque vector
    (*torque)[0] = torques._entity[0][0];
    (*torque)[1] = torques._entity[1][0];
    (*torque)[2] = torques._entity[2][0];
    (*torque)[3] = torques._entity[3][0];   
}