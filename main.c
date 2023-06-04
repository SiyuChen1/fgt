#include <stdio.h>
#include <limits.h>
#include <math.h>

/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
int nchoosek(int n , int k)
{
    int i , n_k = n - k , nchsk = 1;
    if (k < n_k)
    {
        k   = n_k;
        n_k = n - k;
    }
    for ( i = 1 ; i <= n_k ; i++)
    {
        nchsk *= (++k);
        nchsk /= i;
    }
    return nchsk;
}

/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/

double ddist(double *x , double *y , int d)
{
    int i;
    register double t , s = 0.0;
    for (i = 0 ; i < d ; i++)
    {
        t  = (x[i] - y[i]);
        s += (t * t);
    }
    return s;
}
/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/

int idmax(double *x , int N)
{
    int i , k = 0;
    double t = -1.0;

    for (i = 0 ; i < N ; i++ )
    {
        if( t < x[i] )
        {
            t = x[i];
            k = i;
        }
    }
    return k;
}

void Kcenter(double *x , int d , int Nx , int K ,
             double *xc , int *indxc , int *indx , double *xboxsz ,
             double *dist_C)
{
    double *x_ind , *x_j;
    register double temp ;
    int i , j , ind , nd , ibase;

    // randomly pick one node as the first center.
    //	srand( (unsigned)time( NULL ) );
    //	ind      = rand() % Nx;

    ind      = 1;
    *indxc++ = ind;
    x_j      = x;
    x_ind    = x + ind*d;

    for (j = 0 ; j < Nx ; x_j += d , j++)
    {
        dist_C[j] = (j==ind) ? 0.0 : ddist(x_j , x_ind , d);
        printf("current index = %d, distance = %f\n", j, dist_C[j]);
        indx[j]   = 0;
    }
    for(i = 1 ; i < K ; i++)
    {
        ind      = idmax(dist_C , Nx);
        *indxc++ = ind;
        x_j      = x;
        x_ind    = x + ind*d;
        for (j = 0 ; j < Nx ; x_j += d, j++)
        {
            temp = (j==ind) ? 0.0 : ddist(x_j , x_ind , d);
            if (temp < dist_C[j])
            {
                dist_C[j] = temp;
                indx[j]   = i;
            }
        }
    }

    for(int i = 0; i < Nx; i++){
        printf("current index = %d, center_id = %d\n", i, indx[i]);
    }

    for (i = 0 ; i < K ; i++)
    {
        xboxsz[i] = 0;
    }
    for (i = 0; i < d*K; i++)
    {
        xc[i] = 0.0;
    }

    for (i = 0 , nd = 0 ; i < Nx ; i++ , nd += d)
    {
        xboxsz[indx[i]]++;
        ibase = indx[i]*d;
        for (j = 0 ; j < d; j++)
        {
            xc[j + ibase ] += x[j + nd];
        }
    }

    for (i = 0 , ibase = 0 ; i < K ; i++ , ibase += d)
    {
        temp = 1.0/xboxsz[i];
        for (j = 0; j < d; j++)
        {
            xc[j + ibase] *= temp;
        }
    }
}

void Compute_C_k(int d , int p ,
                 double *C_k ,
                 int *heads , int *cinds)
{
    int i , k , t , j , tail , head;

    for (i = 0; i < d; i++)
    {
        heads[i] = 0;
    }

    heads[d] = INT_MAX;
    cinds[0] = 0;
    C_k[0]   = 1.0;

    for (k=1 , t=1, tail = 1 ; k < p ; k++ , tail=t)
    {
        for (i = 0; i < d; i++)
        {
            head     = heads[i];
            heads[i] = t;

            for ( j = head ; j < tail ; j++ , t++)
            {
                cinds[t] = (j < heads[i+1]) ? cinds[j] + 1 : 1;
                C_k[t]   = 2.0 * C_k[j];
                C_k[t]  /= (double) cinds[t];
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------*/
void Compute_A_k(double *x , double *w , double *xc, double *C_k , double sigma , int d , int Nx , int p , int K , int pd ,
                 double *A_k ,
                 int *indx  , double *dx , double *prods , int *heads )
{

    int n , i , k , t , tail , j , head , ind;
    int nbase , ix2c , ix2cbase;
    register double sum , ctesigma = 1.0/(sigma) , temp , temp1;

    for (i = 0; i < pd*K; i++)
    {
        A_k[i] = 0.0;
    }

    for (n = 0 ; n < Nx ; n++)
    {
        printf("n = %d\n", n);
        nbase    = n*d;
        ix2c     = indx[n];
        ix2cbase = ix2c*d;
        ind      = ix2c*pd;
        temp     = w[n];
        sum      = 0.0;
        for (i = 0 ; i < d ; i++)
        {
            dx[i]    = (x[i + nbase] - xc[i + ix2cbase])*ctesigma;
            sum     += dx[i] * dx[i];
            heads[i] = 0;
        }
        prods[0] = exp(-sum);

        printf("before: \n");
        for (i = 0 ; i < pd ; i++)
        {
            printf("prods[%d] = %f\t", i, prods[i]);
        }

        for (k = 1 , t = 1 , tail = 1 ; k < p ; k++ , tail = t)
        {
            for (i = 0 ; i < d; i++)
            {
                head     = heads[i];
                heads[i] = t;
                temp1    = dx[i];

                for ( j = head; j < tail ; j++, t++)
                {
                    printf("j = %d, tail = %d, t = %d\n", j, tail, t);
                    prods[t] = temp1 * prods[j];
                }
            }
        }

        printf("after: \n");
        for (i = 0 ; i < pd ; i++)
        {
            printf("prods[%d] = %f\t", i, prods[i]);
        }
        printf("\n");

        for (i = 0 ; i < pd ; i++)
        {
            A_k[i + ind] += temp*prods[i];
        }


        for (k = 0 ; k < K ; k++)
        {
            ind  = k*pd;
            for (i = 0 ; i < pd ; i++)
            {
                printf("A_k = %f\n", A_k[i + ind]);
            }
        }
    }

    for (k = 0 ; k < K ; k++)
    {
        ind  = k*pd;
        for (i = 0 ; i < pd ; i++)
        {
            printf("A_k = %f\n", A_k[i + ind]);
        }
    }
    printf("refresh\n");

    for (k = 0 ; k < K ; k++)
    {
        ind  = k*pd;
        for (i = 0 ; i < pd ; i++)
        {
            A_k[i + ind] *= C_k[i];
        }
    }
}


void fgt_predict(double *y , double *xc , double *A_k  , int Ny, double sigma , int K , int p, double e , int d , int pd ,
                 double *v ,
                 double *dy , double *prods , int *heads)
{
    int i , j , m , k , t , tail , mbase , kn , xbase , head , ind;
    double sum2 , ctesigma = 1.0/(sigma) , temp , temp1;

    for (m=0 ; m < Ny ; m++)
    {
        temp    = 0.0;
        mbase   = m*d;
        for (kn = 0 ; kn < K ; kn++)
        {
            xbase = kn*d;
            ind   = kn*pd;
            sum2  = 0.0;
            for (i = 0 ; i < d ; i++)
            {
                dy[i]    = (y[i + mbase] - xc[i + xbase])*ctesigma;
                sum2    += dy[i] * dy[i];
                heads[i] = 0;
            }
            if (sum2 > e) continue; //skip to next kn

            prods[0] = exp(-sum2);


            for (k=1, t=1, tail=1 ; k < p ; k++ , tail=t)
            {
                for (i = 0 ; i < d; i++)
                {
                    head     = heads[i];
                    heads[i] = t;
                    temp1    = dy[i];
                    for (j = head ; j < tail ; j++ , t++)
                    {
                        prods[t] = temp1 * prods[j];
                    }
                }
            }

            printf("valid");


            for (i = 0 ; i < pd ; i++)
            {
                printf("%f\n", prods[i]);
                temp += A_k[i + ind]*prods[i];
            }
            printf("sum = %f\n", temp);
        }
        v[m] = temp;
    }
}

int main() {
    printf("Hello, World!\n");
    const int D = 2;
    const int Nx = 5;
    const int K = 2;
    double xc[Nx];
    int indxc[Nx], indx[Nx];
    double xboxsz[Nx * D];
    double dist_C[Nx * D];
    double x[Nx * D];
    x[0] = 10;
    x[1] = 12;
    x[2] = 12;
    x[3] = 14;
    x[4] = 13;
    x[5] = 29;
    x[6] = 14;
    x[7] = 32;
    x[8] = 15;
    x[9] = 11;
    Kcenter(x, D, Nx, K, xc, indxc, indx, xboxsz, dist_C);

    for (int k = 0 ; k < K ; k++)
    {
        for (int d = 0; d < D; d++){
            printf("xc = %f\n", xc[k*D + d]);
        }
    }

    double C_k[100];
    const int P = 6;
    const int R_PD = nchoosek(P + D - 1 , D);
    int heads[100], cinds[100];
    Compute_C_k(D, P, C_k, heads, cinds);

    for (int j = 0 ; j < R_PD ; j++)
    {
        printf("C_k = %f\n", C_k[j]);
    }

    double prods[100], w[Nx], A_k[100], dx[D];
    w[0] = 10;
    w[1] = 10;
    w[2] = 91;
    w[3] = 51;
    w[4] = 21;
    double sigma = 2;
    double pd = nchoosek(P + D - 1 , D);
    Compute_A_k(x, w, xc, C_k, sigma, D, Nx, P, K, pd, A_k, indx, dx, prods, heads);

    for (int k = 0 ; k < K ; k++)
    {
        int ind  = k*pd;
        for (int i = 0 ; i < pd ; i++)
        {
            printf("A_k = %f\n", A_k[i + ind]);
        }
    }

    int Ny = 3;
    double y[Ny * D];
    double v[Ny];
    double dy[D];
    double e = 10;
    y[0] = 14;
    y[1] = 23;
    y[2] = 17;
    y[3] = 32;
    y[4] = 19;
    y[5] = 45;
    fgt_predict(y, xc, A_k, Ny, sigma, K, P, e, D, R_PD, v, dy, prods, heads);
    for (int i = 0 ; i < Ny ; i++)
    {
        printf("v = %f\n", v[i]);
    }
    return 0;
}
