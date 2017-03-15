#ifndef PEXIT_AWGN_HPP
#define PEXIT_AWGN_HPP

#include <cmath>
#include "ProtographClass.hpp"

#define PEXIT_AWGN_TOLERANCE 0.000001

double PEXIT_Threshold_AWGN(const Protograph &P, const double &EbN0_linear_min, const double &EbN0_linear_max, const unsigned int &max_its);
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its, unsigned int &its_taken, double &MI_reached);
bool PEXIT_AWGN(const Protograph &P, const double &EbN0_linear, const unsigned int &max_its);
void PEXIT_AWGN_Initialize(const Protograph &P, const Eigen::VectorXi<double> EbN0_linear, const double R, Eigen::VectorXi<double> &sigma_ch_squared, Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV, Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC, Eigen::VectorXi<double> &IAPP);
void PEXIT_AWGN_Initialize(const Protograph &P, const double EbN0_linear, const double R, Eigen::VectorXi<double> &sigma_ch_squared, Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV, Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC, Eigen::VectorXi<double> &IAPP);
void PEXIT_AWGN_IEV(const Protograph &P, const Eigen::VectorXi<double> &sigma_ch_squared, const Eigen::MatrixXi<double> &IAV, Eigen::MatrixXi<double> &IEV);
void PEXIT_AWGN_IEC(const Protograph &P, const Eigen::MatrixXi<double> &IAC, Eigen::MatrixXi<double> &IEC);
bool PEXIT_AWGN_IAPP(const Protograph &P, const Eigen::VectorXi<double> &sigma_ch_squared, const Eigen::MatrixXi<double> &IAV, Eigen::VectorXi<double> &IAPP, double &IAPP_avg, double &IAPP_min);

#endif