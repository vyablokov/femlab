#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double** A;
double** An;
const double botL = 8.0;
const double leftL = 6.0;
const double rightL = 3.0;
const double topL = 5.0;
const double radius = 3.0;
int Nx, Ny, Nt;
double hx, hy, ht;
double a=0.2;
double leftTemp = 200;
double rightTemp = 100;
double botTemp = 50;
double topTemp = 100;
double edgeTemp = 100;
int bndJ, bndK;
int stepNo;
double k1=1, k2=1, k12=1;

int getEdgeK(int j)
{
    if(j == Nx)
        return bndK;
    if(j > bndJ)
    {
        int edgeK = ceil((sqrt(radius*radius - ((double)j*hx - (double)bndJ*hx)*((double)j*hx - (double)bndJ*hx)) + (double)bndK*hy) / hy);
        return edgeK;
    }
    else
        return Ny;
}

int getEdgeJ(int k)
{
    if(k == Ny)
        return bndJ;
    if(k > bndK)
    {
        int edgeJ = ceil((sqrt(radius*radius - ((double)k*hy - (double)bndK*hy)*((double)k*hy - (double)bndK*hy)) + (double)bndJ*hx) / hx);
        return edgeJ;
    }
    else
        return Nx;
}

double getXEdgeCoeff(int k) {
    int lastXStep = getEdgeJ(k) - 1;
    double edgeCoeffX = sqrt(radius*radius - (k - bndK)*(k - bndK)*hy*hy)/hx - lastXStep + bndJ;
    return edgeCoeffX;
}

double getYEdgeCoeff(int j) {
    int lastYStep = getEdgeK(j) - 1;
    double edgeCoeffY = sqrt(radius*radius - (j - bndJ)*(j - bndJ)*hx*hx)/hy - lastYStep + bndK;
    return edgeCoeffY;
}


double getNext(int j, int k, double** Tn, double** T)
{
    if (j == 0 && (k == Ny / 2) ) {
        return T[j+1][k]/(1 + a*hx);
    }
    if (j == 0 || k == 0)
        return T[j][k];
    double d2T_dx2, d2T_dy2;
    if (j == getEdgeJ(k) - 1 && k > bndK)
    {
        double edgeCoeffX = getXEdgeCoeff(k);
        d2T_dx2 = 2*(T[j+1][k] - (edgeCoeffX+1)*T[j][k] + edgeCoeffX*T[j-1][k]) / (edgeCoeffX*(edgeCoeffX+1)*hx*hx);
    }
    else
        d2T_dx2 = (T[j+1][k] - 2*T[j][k] + T[j-1][k]) / (hx*hx);

    if (j > bndJ &&k == getEdgeK(j+1) - 1)
    {
        double edgeCoeffY = getYEdgeCoeff(j);
        d2T_dy2 = 2*(T[j][k+1] - (edgeCoeffY+1)*T[j][k] + edgeCoeffY*T[j][k-1]) / (edgeCoeffY*(edgeCoeffY+1)*hy*hy);
    }
    else
        d2T_dy2 = (T[j][k+1] - 2*T[j][k] + T[j][k-1]) / (hy*hy);

    return a*ht*(d2T_dx2 + d2T_dy2) + T[j][k];
}

void makeStep(double** newMx, double** oldMx)
{
    for (int j = 0; j < Nx; j++){
        for (int k = 0; k < getEdgeK(j+1); k++){
            newMx[j][k] = getNext(j, k, newMx, oldMx);
        }
    }
}

void printMx(double** mat)
{
    FILE* s = fopen("mx.txt", "w");

    for (int k = Ny; k >= 0; k--){
        for (int j = 0; j <= Nx; j++)
            fprintf(s, "%.4f\t", mat[j][k]);
        fputc('\n', s);
    }
    putchar('\n');
    fclose(s);
}

void printPoints(double** mat, int num)
{
    FILE* s = fopen("points.txt", "w");
    double x, y;
    double edgeX, edgeY;
    int j,k;
    int edgeJ, edgeK;
    for (j=0; j <= Nx; j++) {
        edgeK = getEdgeK(j);
        for(k = 0; k < edgeK; k++) {
            edgeJ = getEdgeJ(k);
            x = hx*j;
            y = hy*k;
            if(j == (edgeJ - 1) && (j > bndJ && k > bndK)) {
                edgeX = hx*(j+getXEdgeCoeff(k));
                fprintf(s, "%.2f %.2f %.2f\n", edgeX, y, mat[j+1][k]);
            }
            if(k == (edgeK - 1) && (j > bndJ && k > bndK)) {
                edgeY = hy*(k+getYEdgeCoeff(j));
                fprintf(s, "%.2f %.2f %.2f\n", x, edgeY, mat[j][k+1]);
            }
            fprintf(s, "%.2f %.2f %.2f\n", x, y, mat[j][k]);
        }
        if (j <= bndJ || j == Nx)
            fprintf(s, "%.2f %.2f %.2f\n", x, edgeK*hx, mat[j][Ny]);
        fprintf(s, "\n");
    }
    fclose(s);
    system("gnuplot \"plot.sc\"");
}

void clearMx(double** A)
{
    for (int j = 0; j < Nx+1; j++)
        for (int k = 0; k < Ny+1; k++)
            A[j][k] = 0;
}

void copyMx(double** dst, double** src)
{
    for (int j = 0; j < Nx+1; j++)
        for (int k = 0; k < Ny+1; k++)
            dst[j][k] = src[j][k];
}

int main(int argc, char *argv[])
{
    if (argc < 4) {
        printf("Usage: ./femlab <X steps> <Y steps> <N steps>\n");
        exit(1);
    }
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    ht = 0.001;
    Nt = atoi(argv[3]);

    hx = botL / Nx;
    hy = leftL / Ny;
    bndJ = (int)Nx - (botL - topL)/hx;
    bndK = (int)Ny - (leftL - rightL)/hy;

    A = new double*[Nx+1];
    for (int i = 0; i < Nx+1; i++)
        A[i] = new double[Ny+1];

    for (int j = 0; j < Nx+1; j++){
        for (int k = 0; k < Ny+1; k++){
            A[j][k] = 0.0;
        }
    }

    An = new double*[Nx+1];
    for (int i = 0; i < Nx+1; i++)
        An[i] = new double[Ny+1];

    for (int j = 0; j < Nx+1; j++){
        for (int k = 0; k < Ny+1; k++){
            An[j][k] = 100.0;
        }
    }

    for (int j =0; j < Nx+1; j++){
        A[j][0] = botTemp;
        An[j][0] = botTemp;
    }
    for (int k = 0; k < Ny+1; k++){
        A[0][k] = leftTemp;
        An[0][k] = leftTemp;
    }

    for (int j = 0; j < Nx; j++) {
        for (int k = getEdgeK(j); k >= getEdgeK(j+1); k--) {
            A[j][k] = topTemp;
            An[j][k] = topTemp;
        }
    }

    for (int k = 0; k < bndK+1; k++){
        A[Nx][k] = rightTemp;
        An[Nx][k] = rightTemp;
    }

    stepNo = 0;
    for (stepNo=0; stepNo < Nt; stepNo++){
        makeStep(An, A);
        copyMx(A, An);
    }
    printMx(A);
    printPoints(A, stepNo);
    printf("Ready.\n");

    return 0;
}
