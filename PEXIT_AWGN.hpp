#ifndef PEXIT_AWGN_HPP
#define PEXIT_AWGN_HPP

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include "ProtographClass.hpp"

/* NB, we need 64 bit capability for this */

#define PEXIT_AWGN_TOLERANCE 0.000001
#define JLUT_MULTIPLIER 1000000
#define JINVLUT_MULTIPLIER 1000000
#define SIGMA_STAR 1.6363
#define I_STAR 0.3646

#define A_J1 -0.0421061
#define B_J1 0.209252
#define C_J1 -0.00640081

#define A_J2 0.00181491
#define B_J2 -0.142675
#define C_J2 -0.0822054
#define D_J2 0.0549608

#define A_S1 1.09542
#define B_S1 0.214217
#define C_S1 2.33727

#define A_S2 0.706692
#define B_S2 0.386013
#define C_S2 -1.75017

double PEXIT_Threshold_AWGN(const Protograph &P, const double &EbN0_linear_min, const double &EbN0_linear_max, const unsigned int &max_its);
double PEXIT_Sweep_AWGN(const Protograph &P, const double &EbN0_linear_min, const double &EbN0_linear_max, const double &step, const unsigned int &max_its);
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its, unsigned int &its_taken, double &MI_reached);
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its);
void PEXIT_AWGN_Initialize(const Protograph &P, const Eigen::VectorXd EbN0_linear, const double R, Eigen::VectorXd &sigma_ch_squared, Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV, Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC, Eigen::VectorXd &IAPP);
void PEXIT_AWGN_Initialize(const Protograph &P, const double EbN0_linear, const double R, Eigen::VectorXd &sigma_ch_squared, Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV, Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC, Eigen::VectorXd &IAPP);
void PEXIT_AWGN_IEV(const Protograph &P, const Eigen::VectorXd &sigma_ch_squared, const Eigen::MatrixXd &IAV, Eigen::MatrixXd &IEV);
void PEXIT_AWGN_IEC(const Protograph &P, const Eigen::MatrixXd &IAC, Eigen::MatrixXd &IEC);
bool PEXIT_AWGN_IAPP(const Protograph &P, const Eigen::VectorXd &sigma_ch_squared, const Eigen::MatrixXd &IAV, Eigen::VectorXd &IAPP, double &IAPP_avg, double &IAPP_min);

double dBtoLin(const double &);
double lintodB(const double &);

void createJLUT(std::vector<double> &);
double getJ(const double &, const std::vector<double> &);
double getJ(const double &);

void createJinvLUT(std::vector<double> &);
double getJinv(const double &, const std::vector<double> &);
double getJinv(const double &);

#endif