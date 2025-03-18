#include <dace/dace.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_multiroots.h>

#include <fmt/printf.h>
#include <gsl/gsl_vector.h>
#include "RK78.h"
#include "AstroCnvRef.h"
#include "constants.h"
#include <chrono>
#include "osctomean.h"
#include "convert.h"
#include "dynamics.h"
#include <eigen/Eigen/Dense>
#include "calcnonlinearindx.h"
#include "uncertainityConvert.h"
 
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

        cst[9] = inputs[0][6]; //params.LU;
        cst[10] = inputs[0][7]; //params.Area;
        cst[11] = inputs[0][8]; // params.Cd;
        cst[12] = inputs[0][9];//params.MU;
        cst[13] = inputs[0][13]; //param.drag

        const int Nnodes = inputs[0][10];
        const double tf = inputs[0][11];
        const double nu = inputs[0][12];
        const double coordSys = inputs[0][14];

   
        typedef DACE::DA state_type; // define state type as DA
        DACE::DA::init(2, 6);        // initialise DACE for X order in N variables
        DACE::DA::setEps(1e-40);
        DACE::AlgebraicVector<state_type> x0(6), x0t(6), fx(6), fxbar(6);
  
        //DACE::AlgebraicMatrix<double> ;
        double statesGuess[6], xarray[6], constantPart[6],dx[6], dx0[6];
        Eigen::MatrixXd Mat(6,6);
        AlgebraicVector<state_type> EigenV(6,0.0);
        AlgebraicVector<state_type> Geqoe_Eigen_DA(6);
        AlgebraicVector<state_type> Geqoesc(6);
        DACE::AlgebraicMatrix<double> STM1(6,6), STM0(6,6);
        DACE::AlgebraicMatrix<state_type> STM(6,6);  

        vector <double> dxs(6), dx0s(6), fxe(6),fx0(6);
    

        ArrayFactory f;

        // Cartesian

        outputs[0] =f.createArray<double>({6});
        outputs[1] =f.createArray<double>({6});
        outputs[2] =f.createArray<double>({6,6});
        outputs[3] =f.createArray<double>({6,6});


        // MEE 
        outputs[4] =f.createArray<double>({6});
        outputs[5] =f.createArray<double>({6});
        outputs[6] =f.createArray<double>({6,6});
        outputs[7] =f.createArray<double>({6,6});


        // Geqoe
        outputs[8] =f.createArray<double>({6});
        outputs[9] =f.createArray<double>({6});
        outputs[10] =f.createArray<double>({6,6});
        outputs[11] =f.createArray<double>({6,6});



        // Ceqoe
        outputs[12] =f.createArray<double>({6});
        outputs[13] =f.createArray<double>({6});
        outputs[14] =f.createArray<double>({6,6});
        outputs[15] =f.createArray<double>({6,6});


        // Kepler
        outputs[16] =f.createArray<double>({6});
        outputs[17] =f.createArray<double>({6});
        outputs[18] =f.createArray<double>({6,6});
        outputs[19] =f.createArray<double>({6,6});



        for (unsigned int j = 0; j < 6; j++){
            statesGuess[j] = inputs[1][j];
        }

        for (unsigned int j = 0; j < 6; j++){
            dx[j] = inputs[2][j];
        }

          for (unsigned int j = 0; j < 6; j++){
             dxs[j] = dx[j];
              dx0s[j] = 0.0; 
        }


    
        
        int Jsum = 0;

        for (unsigned int i = 0; i < 6; i++)
        {
            if (i <3)
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
            else
            {
                x0[i] = statesGuess[i] + DACE::DA(i + 1);
            }
        }

        // propagate x0, and output the result please.
        fx = RK78(6, x0, 0.0, tf, cst, DynamicsCart);

        STM=  getSTM(tf, fx);

        fxe = fx.eval(dxs);
        fx0 = fx.eval(dx0s);
       
        for (unsigned int i = 0; i < 6; i++){
            outputs[0][i] = fxe[i];
            outputs[1][i] = fx0[i]; 
            for (unsigned int j = 0; j < 6; j++)
            {
                STM0.at(i,j) = STM.at(i,j).eval(dx0s);
                STM1.at(i,j) = STM.at(i,j).eval(dxs);

                outputs[2][i][j] = STM0.at(i,j);
                outputs[3][i][j] = STM1.at(i,j);


            }
        }

        // MEE

        //convert to mee
        x0t = eci2mee(x0, cst[5]);
        Geqoesc = convertU(tf, x0t);

        //propagate
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsMEE);
        //fx = eci2mee(fx, cst[5]);

        STM=  getSTM(tf, fx);

        fxe = fx.eval(dxs);
        fx0 = fx.eval(dx0s);

        for (unsigned int i = 0; i < 6; i++){
            outputs[4][i] = fxe[i];
            outputs[5][i] = fx0[i];
            for (unsigned int j = 0; j < 6; j++)
            {
                STM0.at(i,j) = STM.at(i,j).eval(dx0s);
                STM1.at(i,j) = STM.at(i,j).eval(dxs);

                outputs[6][i][j] = STM0.at(i,j);
                outputs[7][i][j] = STM1.at(i,j);


            }
        }


        // CEQoe
        x0t     = RV2CEq(x0, cst);
        Geqoesc = convertU(tf, x0t);

        //propagate
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsCeq);
        //fx = RV2CEq(fx, cst);

        STM=  getSTM(tf, fx);

        fxe = fx.eval(dxs);
        fx0 = fx.eval(dx0s);

        for (unsigned int i = 0; i < 6; i++){
            outputs[8][i] = fxe[i];
            outputs[9][i] = fx0[i];
            for (unsigned int j = 0; j < 6; j++)
            {
                STM0.at(i,j) = STM.at(i,j).eval(dx0s);
                STM1.at(i,j) = STM.at(i,j).eval(dxs);

                outputs[10][i][j] = STM0.at(i,j);
                outputs[11][i][j] = STM1.at(i,j);


            }
        }

        // GEQoe
        x0t     = RV2GEq(x0, cst);
        Geqoesc = convertU(tf, x0t);

        //propagate

        
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsEq);
        //fx = RV2GEq(fx, cst);

        STM=  getSTM(tf, fx);

        fxe = fx.eval(dxs);
        fx0 = fx.eval(dx0s);

        for (unsigned int i = 0; i < 6; i++){
            outputs[12][i] = fxe[i];
            outputs[13][i] = fx0[i];
            for (unsigned int j = 0; j < 6; j++)
            {
                STM0.at(i,j) = STM.at(i,j).eval(dx0s);
                STM1.at(i,j) = STM.at(i,j).eval(dxs);

                outputs[14][i][j] = STM0.at(i,j);
                outputs[15][i][j] = STM1.at(i,j);


            }
        }

        // Kepler
        x0t     = eci2Kep(x0,cst[5]);
        Geqoesc = convertU(tf, x0t);

//         fx = RK78(6, Kep2eci(Geqoesc, cst[5]), 0.0, tf, cst, DynamicsCart);
//         fx = eci2Kep(fx, cst[5]);

    
//            //propagate
        fx = RK78(6, Geqoesc, 0.0, tf, cst,DynamicsKep);
//        // fx = eci2Kep(fx, cst[5]);
//         
        STM=  getSTM(tf, fx);

        fxe = fx.eval(dxs);
        fx0 = fx.eval(dx0s);

        for (unsigned int i = 0; i < 6; i++){
            outputs[16][i] = fxe[i];
            outputs[17][i] = fx0[i];
            for (unsigned int j = 0; j < 6; j++)
            {
                STM0.at(i,j) = STM.at(i,j).eval(dx0s);
                STM1.at(i,j) = STM.at(i,j).eval(dxs);

                outputs[18][i][j] = STM0.at(i,j);
                outputs[19][i][j] = STM1.at(i,j);


            }
        }






        








    }
   
    

    };
