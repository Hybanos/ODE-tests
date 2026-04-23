#include "odes.hpp"

// reference DOP853 implementation, ugly code don't look

ref_ftype haha;
std::ofstream fout;
int ndims;
bool _save_results;

void dop_ref_solout(int *nr, double *xold, double *X, double *Y, 
                    int *n, double *con, int *icomp, int *nrd, 
                    int *rpar, int *ipar, int *irtrn, double *xout) {

    if (!_save_results) return;

    int bodies = ndims / 7;
    vecarray x = vecarray((vec3 *) Y, bodies);
    vecarray v = vecarray((vec3 *) Y + bodies, bodies);
    array m = array(Y + bodies * 6, bodies);

    fout << std::setprecision(10)
      << bodies << ";" << *nr << ";" << *X << ";";
    for (int i = 0; i < bodies; i++) {
        fout << m[i] << ";" 
          << x[i].x << ";" << x[i].y << ";" << x[i].z << ";"
          << v[i].x << ";" << v[i].y << ";" << v[i].z << ";";
    }

    double gamma = 1.0;
    double K = 0.0;
    double U = 0.0;
    for (int i = 0; i < bodies; i++) {
        K += v[i].norm() * v[i].norm() * m[i] / 2.0;
        for (int j = i+1; j < m.size(); j++) {
            U += -gamma * m[i] * m[j] / (x[j] - x[i]).norm();
        }
    }

    vec3 A = vec3{0, 0, 0}; 
    for (int i = 0; i < bodies; i++) {
        A = A + v[i].prod(x[i]);
    }

    fout << A.norm() << ";" << K << ";" << U << ";" << K+U
      << std::endl;
}

void dop_ref_f_wrapper(int *n, double *t, double * X, double * F, double *idk, int *idk2) {
    haha(n, t, X, F, idk, idk2);
}

void DOP853_ref::step() {
    int n = Y.extent(0);
    ndims = n;
    _save_results = this->save_results;
    haha = [&](int *n, double *t, double * X, double * F, double *idk, int *idk2){
        array x(X, *n);
        array f(F, *n);
        this->f(this->t, x, f);
    };
    double x = 0.0;
    double *y = Y.data_handle();
    double xend = 10.0;
    double rtol = R_TOL;
    double atol = A_TOL;
    int itol = 0;
    int iout = 1;
    int lwork = 11 * n + 21;
    double *work = new double[lwork];
    int liwork = 21;
    int *iwork = new int[liwork];
    double rpar = 0;
    int ipar = 0;
    int idid = 0;

    for (int i = 0; i < lwork; i++) work[i] = 0.0;
    for (int i = 0; i < liwork; i++) iwork[i] = 0;

    // std::cout << (void *) dop_ref_f_wrapper << std::endl;

    fout.open("DOP853-ref.txt", std::ios::app);
    dop853(&n, (void *)dop_ref_f_wrapper, &x, Y.data_handle(), 
        &xend, &rtol, &atol, &itol, (void *)dop_ref_solout, 
        &iout, work, &lwork, iwork, &liwork, 
        &rpar, &ipar, &idid);
    fout.close();
    t = 10;
    // std::cout << t << std::endl;
    f_evals = iwork[16];
    steps = iwork[18];
}