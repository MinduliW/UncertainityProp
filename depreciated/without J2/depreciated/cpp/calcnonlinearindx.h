#ifndef calcnonlinearindx_H
#define calcnonlinearindx_H


template <typename T, typename U>
U nonlinearindex(U t, DACE::AlgebraicVector<T> fx, DACE::AlgebraicVector<T> fxbar)
{
    DACE::AlgebraicMatrix<U> Jacobian(6, 6), Jacobianbar(6,6),diff(6,6), xval;
    DACE::AlgebraicVector<U> colsum(6),colsumJbar(6);
    U nus;

   // cout << fx.linear() << endl;
    
    
    Jacobian = fx.linear();
    Jacobianbar = fxbar.linear();

    diff = Jacobian - Jacobianbar;

    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 0; j < 6; j++)
        {
            colsum[i] += DACE::cons(diff.at(i,j));
            colsumJbar[i] += DACE::cons(Jacobianbar.at(i,j));
            
        }
    
    }

    nus = *max_element(colsum.begin(), colsum.end())/  *max_element(colsumJbar.begin(), colsumJbar.end());
   // cout << nus << endl;

    //nus = max(colsum)/max(colsumJbar);


    return nus;

}



// 
// template <typename T, typename U>
// DACE::AlgebraicVector<T> nonlinearindex(U t, DACE::AlgebraicVector<T> fx)
// {
//     vector<U> zerostate(6);
//     DACE::AlgebraicMatrix<T> Jacobian(6, 6), xval;
//     DACE::AlgebraicMatrix<U>  Jbar(6, 6);
//     vector<T> zerostateDA(6);
//     DACE::AlgebraicVector<T> nonlinindx(6);
//     T b= 0; 
// 
// 
//     // calculate the jacobian and jbar norm
//     U normJbar = 0;
//     for (unsigned int i = 0; i < 6; i++)
//     {
//         for (unsigned int j = 1; j < 7; j++)
//         {
//             Jacobian.at(i, j - 1) = fx[i].deriv(j);
//             Jbar.at(i, j - 1) = fx[i].deriv(j).eval(zerostate);
//             normJbar += pow(Jbar.at(i, j - 1),2);
//         }
//     }
// 
// 
//     // analyse the dependency on each variable individually
//     for (int l = 0; l < 6; l++)
//     {
//         zerostateDA = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
//         zerostateDA[l] = DACE::DA(l + 1);
//         b = 0;
//         for (unsigned int i = 0; i < 6; i++)
//         {
//             for (unsigned int j = 1; j < 7; j++)
//             {
//                 b += pow(Jacobian.at(i, j - 1).eval(zerostateDA).deriv(l + 1),2);
//             }
//         }
// 
// 
//            // cout <<b << endl;
//         nonlinindx[l] =   (sqrt(b)/sqrt(normJbar));
//       
//     }
// 
//     return nonlinindx;
// 
// }


#endif