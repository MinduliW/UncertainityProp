#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>
#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "constants.h"
#include <chrono>
#include "dynamics.h"
#include "osctomean.h"
#include "convert.h"


#include "mex.hpp"
#include "mexAdapter.hpp"


using namespace std;
using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
    public:
    void operator()(ArgumentList outputs, ArgumentList inputs) override {

        double cst[15];

        cst[3] = inputs[0][0]; //Tmax
        cst[4] = inputs[0][1]; //m0
        cst[5] = inputs[0][2]; //mu
        cst[6] = inputs[0][3]; // J2
        cst[7] = inputs[0][4]; //Re


        const int method = inputs[0][5];


        if (method < 2)
        {
            cst[8] = 1.0;
            cout << "Method 1 running" << endl;
        }
        else
        {
            cst[8] = 2.0;
            cout << "Method 2 running " << endl;
        }

        cst[9] = inputs[0][6]; //params.LU;
        cst[10] = inputs[0][7]; //params.Area;
        cst[11] = inputs[0][8]; // params.Cd;
        cst[12] = inputs[0][9];//params.MU;
        cst[13] = inputs[0][13]; //param.drag

        const int Nnodes = inputs[0][10];
        const double dt = inputs[0][11];
        const double nu = inputs[0][12];
        const double coordSys = inputs[0][14];

        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 9);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), fx(6), dv(3), dx(3),fxnew(6),fxendRV(6),fxendkep(6);
        DACE::AlgebraicMatrix<state_type> Jacobian(6, 9), xval(Nnodes, 9);
        DACE::DA derivative(1), result(1), b;
        DACE::AlgebraicMatrix<double> STMgamma(6, 9), STMxftransform(6,9), Jbar(6, 9);
        double statesGuess[Nnodes+ 1][9], xarray[9], constantPart[6], B, normJbar, time, maxval, eclipse;
        vector<double> zerostate(9);
        vector<state_type> zerostateDA(9);
        double dxCeq[6]; 
 

        for (unsigned int i = 0; i < Nnodes + 1; i++){
            for (unsigned int j = 0; j < 9; j++){
                statesGuess[i][j] = inputs[1][i][j];
            }

        }


        // organize outputs
        ArrayFactory f;
        outputs[0] =f.createArray<double>({6*static_cast<unsigned long>(Nnodes),10});
        outputs[1] = f.createArray<double>({static_cast<unsigned long>(Nnodes), 9});
        outputs[2] = f.createArray<double>({6, 7});
        outputs[3] = f.createArray<double>({static_cast<unsigned long>(Nnodes), 6});

        int Jsum = 0; 
        for (int k = 0; k < Nnodes; k++)

        {
            // update time.
            time += dt;

            for (unsigned int i = 0; i < 6; i++)
            {
                x0[i] = statesGuess[k][i] + DACE::DA(i + 1);

                if (i < 3)
                    cst[i] = statesGuess[k][i + 6];
            }


            if (coordSys < 2){
                x0 = GEq2RV(x0, cst);
            }
            else{
                // convert from CEQ to cartesian
                x0 = CEq2RV(x0, cst);
            }


            if (cst[8] < 2)
            {
                fx = RK78(6, x0, 0.0, dt, cst, DynamicsCart);
            }
            else
            {
                for (unsigned int i = 0; i < 3; i++)
                {
                    dv[i] =  (cst[i] + DACE::DA(i + 7)); //
                    x0[i + 3] += dv[i];
                }

                fx = RK78(6, x0, 0.0, dt, cst, DynamicsCart);
            }

            if (coordSys < 2){

                fx = RV2GEq(fx, cst);
            }
            else{
                // convert from CEQ to cartesian

                fx = RV2CEq(fx, cst);
            }

           
            STMgamma = fx.linear();


            if (k == Nnodes-1){

                // do cons(fx)+DA
                for (unsigned int i = 0; i < 6; i++)
                {
                    fxnew[i] = cons(fx[i])+ DACE::DA(i + 1);
                }

                // convert to mean keplerian
                if (coordSys < 2){
                    fxendRV  = GEq2RV(fxnew,cst);
                }
                else{
                    fxendRV  = CEq2RV(fxnew,cst);
                }


                
                fxendRV = oscCart2MeanCart(fxendRV,cst[5], cst[6], cst[7]);
                fxendkep = cart2kep(fxendRV,cst[5]);


                // output the maps.
                STMxftransform = fxendkep.linear();

                for (unsigned int j = 0; j < 6; j++) // Prints row of x
                {

                    for (unsigned int i = 0; i < 6; i++)
                    {
                        outputs[2][j][i] = STMxftransform.at(j, i);
                    }
                    outputs[2][j][6] = cons(fxendkep[j]);
                }
            }

            for (unsigned int i = 0; i < 6; i++)
                constantPart[i] = DACE::cons(fx[i]);

            // output stm and constant term
            for (unsigned int j = 0; j < 6; j++) // Prints row of x
            {
                for (unsigned int i = 0; i < 9; i++)
                    outputs[0][Jsum +j][i]= STMgamma.at(j, i);

                outputs[0][Jsum +j][9] = constantPart[j];
            }

            Jsum += 6;


            // calculate the jacobian and jbar norm
            normJbar = 0;
            for (unsigned int i = 0; i < 6; i++)
            {
                for (unsigned int j = 1; j < 10; j++)
                {
                    Jacobian.at(i, j - 1) = fx[i].deriv(j);
                    Jbar.at(i, j - 1) = fx[i].deriv(j).eval(zerostate);
                    normJbar += pow(Jbar.at(i, j - 1),2);
                }
            }

             // convert some coordinates to Eqoe to get the deviation needed. 
             if (coordSys >= 2.0) {

                 fx = CEq2RV(fx, cst);
                 fx = RV2GEq(fx, cst); 

                 for (unsigned int i = 0; i < 6; i++)
                {
                    fxnew[i] = cons(fx[i])+ DACE::DA(i + 1);
                }

                 fx = GEq2RV(fxnew, cst);
                 fx = RV2CEq(fx, cst); 

                 //cout << fx[5] << endl;

                
                 for (unsigned int i = 0; i < 6; i++)
                 {
                      dxCeq[i]  = 0.0; 
                     for (unsigned int j = 1; j < 7; j++)
                     {
                         dxCeq[i] += abs(cons(fx[i].deriv(j)));
                        // cout <<  cons(fx[i].deriv(j)) << endl;
                      
                     }
                 }

                 //cout << dxCeq<< endl;


             }

            // analyse the dependency on each variable individually
            for (int l = 0; l < 9; l++)
            {
                zerostateDA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                zerostateDA[l] = DACE::DA(l + 1);
                b = 0;
                for (unsigned int i = 0; i < 6; i++)
                {
                    for (unsigned int j = 1; j < 10; j++)
                    {
                        b += pow(Jacobian.at(i, j - 1).eval(zerostateDA).deriv(l + 1),2);
                    }
                }

                xval.at(k, l) =  nu *sqrt(normJbar) / sqrt(b);
                outputs[1][k][l] = DACE::cons(xval.at(k, l));


                if (l < 6){
                // output the nonlinearity index, with unitary variation in Geqoe elements. 
                if (coordSys < 2) 
                {
                    outputs[3][k][l] =  DACE::cons(sqrt(b)/sqrt(normJbar));
                }
                else
                {
                    outputs[3][k][l] =  DACE::cons(sqrt(b)/sqrt(normJbar))*dxCeq[l];

                }
                }

            }

        }
    

    }
};