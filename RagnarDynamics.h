#ifndef RagnarDynamics_h
#define RagnarDynamics_h

#include <Arduino.h>

void ragnarMassCentrGrav(
    double theta[4], double dtheta[4], float (*params)[4][8], 
    double parammass[6], double sc[4][6], double (*mass_matrix)[7][7], 
    double (*coriolis_vector)[7], double (*gravity_vector)[7], 
    double sct[8], double g);
void ragnarMassCentrGravf(
    float theta[4], float dtheta[4], float (*params)[4][8], 
    float parammass[6], float sc[4][6], float (*mass_matrix)[7][7], 
    float (*coriolis_vector)[7], float (*gravity_vector)[7], 
    float sct[8], float g);
void ragnarMassvCentrGrav(
    double theta[4], double dtheta[4], double ddq[7], float params[4][8], 
    double parammass[6], double sc[4][6], double (*mass_vector)[7], 
    double (*coriolis_vector)[7], double (*gravity_vector)[7], double sct[8],
    double g);
void ragnarMassvCentrGravf(
    float theta[4], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6], float (*mass_vector)[7], 
    float (*coriolis_vector)[7], float (*gravity_vector)[7], float sct[8],
    float g);
void ragnarTorques(
    double q[7], double dtheta[4], double ddq[7], float params[4][8], 
    double parammass[6], double sc[4][6],double sct[8], double g, 
    double (*torque)[4]);
void ragnarTorquesf(
    float q[7], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6],float sct[8], float g, 
    float (*torque)[4]);
void ragnarTorquesfimp(
    float q[7], float dtheta[4], float ddq[7], float params[4][8], 
    float parammass[6], float sc[4][6],float sct[8], float g, 
    float (*torque)[4], float impednce[3]);
bool ragnarLegMassCentrGrav(
    float dtheta[4], float detazeta[4][2], float params[4][8], 
    float parammass[6], float sc[4][6], float sct[8],float g,
    float scez[4][4], int leg_num,float (*mass_matrix)[3][3],
    float (*coriolis_vector)[3], float (*gravity_vector)[3]);
void ragnarTorquesNf(
    // This is to get the torques with the task space formulation 
    float q[7], float passivef[4][2], float dtheta[4], float ddx[3], 
    float parameter[4][8], float parammass[6], float sc[4][6],float sct[8], 
    float scez[4][4],float g, float (*torque)[4]);
void computedx(
    // This is to get the dx 
    float q[7], float passivef[4][2], float dtheta[4], float parameter[4][8], 
    float sc[4][6],float sct[8], float (*dx)[3]);
#endif