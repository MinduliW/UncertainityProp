#ifndef calcnonlinearindx_H
#define calcnonlinearindx_H



template <typename T, typename U>
DACE::AlgebraicMatrix<T> getSTM(U t, DACE::AlgebraicVector<T> fx)
{
   
    DACE::AlgebraicMatrix<T> Jacobian(6, 6); 


    // calculate the jacobian and jbar norm
    U normJbar = 0;
    for (unsigned int i = 0; i < 6; i++)
    {
        for (unsigned int j = 1; j < 7; j++)
        {
            Jacobian.at(i, j - 1) = fx[i].deriv(j);
        }
    }


    return Jacobian;

}


#endif