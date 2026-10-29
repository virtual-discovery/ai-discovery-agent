/* Classical Heisenberg Monte Carlo for a dilute two-species magnet (construction reference).
 *
 * H = - sum_{i<j} J_ij e_i . e_j   (each pair once, unit vectors, J in meV, J > 0 ferromagnetic)
 * Heat-bath + over-relaxation sweeps.
 *
 * Usage: mc_dilute graph.bin seed n_eq n_meas n_or T1 [T2 ...]
 * graph.bin (little endian): int32 N, int32 species[N] (0 or 1), int32 offs[N+1], int32 nbr[M], double J[M],
 *            then optionally double pos[3N] (Angstrom) and double box[3] (orthorhombic periodic box, Angstrom)
 * Output, one line per temperature (temperatures are run in the given order, each continuing from
 * the previous configuration, so give them in decreasing order):
 *   T  e_per_spin  then for each order parameter o in {A+B, A-B, A, B}: <m^2> <m^4>, then <|m|> for the four,
 *   then (if positions were given) <|m(k_min)|^2> of the total magnetization, averaged over the three smallest
 *   nonzero wave vectors k = 2 pi / box_d along x, y, z (same normalisation as <m^2>)
 * where m is the order-parameter vector divided by the number of spins of the species involved.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t s[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next(void) {
    const uint64_t result = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
}
static inline double uni(void) { return (next() >> 11) * 0x1.0p-53; }
static void seed_rng(uint64_t x) {
    for (int i = 0; i < 4; i++) { x += 0x9e3779b97f4a7c15ULL; uint64_t z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL; z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; s[i] = z ^ (z >> 31); }
}

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage\n"); return 1; }
    FILE *f = fopen(argv[1], "rb");
    int N; if (fread(&N, 4, 1, f) != 1) return 2;
    int *sp = malloc(4 * N), *offs = malloc(4 * (N + 1));
    if (fread(sp, 4, N, f) != (size_t)N || fread(offs, 4, N + 1, f) != (size_t)(N + 1)) return 2;
    int M = offs[N];
    int *nb = malloc(4 * (M > 0 ? M : 1)); double *J = malloc(8 * (M > 0 ? M : 1));
    if (M > 0 && (fread(nb, 4, M, f) != (size_t)M || fread(J, 8, M, f) != (size_t)M)) return 2;
    double *pos = malloc(sizeof(double) * 3 * (N > 0 ? N : 1)), box[3];
    int haspos = fread(pos, 8, 3 * N, f) == (size_t)(3 * N) && fread(box, 8, 3, f) == 3;
    fclose(f);
    double *cs = malloc(sizeof(double) * 6 * (N > 0 ? N : 1));
    if (haspos) for (int i = 0; i < N; i++) for (int d = 0; d < 3; d++) {
        double ph = 2 * M_PI * pos[3*i+d] / box[d]; cs[6*i+2*d] = cos(ph); cs[6*i+2*d+1] = sin(ph); }
    seed_rng(strtoull(argv[2], 0, 10));
    long neq = atol(argv[3]), nmeas = atol(argv[4]); int nor = atoi(argv[5]);
    int nA = 0, nB = 0; for (int i = 0; i < N; i++) { if (sp[i]) nB++; else nA++; }
    double *e = malloc(sizeof(double) * 3 * N);
    for (int i = 0; i < N; i++) { /* start ferromagnetic along z */ e[3*i] = 0; e[3*i+1] = 0; e[3*i+2] = 1; }
    const double kB = 0.08617333262;
    for (int it = 6; it < argc; it++) {
        double T = atof(argv[it]), beta = 1.0 / (kB * T);
        double acc[12] = {0}, eacc = 0, skacc = 0; long nm = 0;
        for (long sweep = 0; sweep < neq + nmeas; sweep++) {
            for (int i = 0; i < N; i++) { /* heat bath */
                double hx = 0, hy = 0, hz = 0;
                for (int k = offs[i]; k < offs[i+1]; k++) { int j = nb[k]; hx += J[k]*e[3*j]; hy += J[k]*e[3*j+1]; hz += J[k]*e[3*j+2]; }
                double h = sqrt(hx*hx + hy*hy + hz*hz), x = beta * h, c;
                if (x < 1e-8) c = 2 * uni() - 1;
                else { double u = uni(); c = 1 + log(u + (1 - u) * exp(-2 * x)) / x; if (c < -1) c = -1; if (c > 1) c = 1; }
                double sn = sqrt(1 - c*c), ph = 2 * M_PI * uni(), a1 = sn * cos(ph), a2 = sn * sin(ph);
                double ux, uy, uz;
                if (x < 1e-8) { ux = 0; uy = 0; uz = 1; } else { ux = hx/h; uy = hy/h; uz = hz/h; }
                /* orthonormal frame around u */
                double vx, vy, vz;
                if (fabs(uz) < 0.9) { vx = -uy; vy = ux; vz = 0; } else { vx = 0; vy = -uz; vz = uy; }
                double vn = sqrt(vx*vx + vy*vy + vz*vz); vx /= vn; vy /= vn; vz /= vn;
                double wx = uy*vz - uz*vy, wy = uz*vx - ux*vz, wz = ux*vy - uy*vx;
                e[3*i] = c*ux + a1*vx + a2*wx; e[3*i+1] = c*uy + a1*vy + a2*wy; e[3*i+2] = c*uz + a1*vz + a2*wz;
            }
            for (int r = 0; r < nor; r++) for (int i = 0; i < N; i++) { /* over-relaxation */
                double hx = 0, hy = 0, hz = 0;
                for (int k = offs[i]; k < offs[i+1]; k++) { int j = nb[k]; hx += J[k]*e[3*j]; hy += J[k]*e[3*j+1]; hz += J[k]*e[3*j+2]; }
                double h2 = hx*hx + hy*hy + hz*hz; if (h2 < 1e-24) continue;
                double p = 2 * (e[3*i]*hx + e[3*i+1]*hy + e[3*i+2]*hz) / h2;
                e[3*i] = p*hx - e[3*i]; e[3*i+1] = p*hy - e[3*i+1]; e[3*i+2] = p*hz - e[3*i+2];
            }
            if (sweep >= neq) {
                double MA[3] = {0}, MB[3] = {0}, en = 0;
                for (int i = 0; i < N; i++) {
                    double *m = sp[i] ? MB : MA; m[0] += e[3*i]; m[1] += e[3*i+1]; m[2] += e[3*i+2];
                    for (int k = offs[i]; k < offs[i+1]; k++) { int j = nb[k]; en -= 0.5 * J[k] * (e[3*i]*e[3*j] + e[3*i+1]*e[3*j+1] + e[3*i+2]*e[3*j+2]); }
                }
                double o[4][3];
                for (int d = 0; d < 3; d++) {
                    o[0][d] = (MA[d] + MB[d]) / N; o[1][d] = (MA[d] - MB[d]) / N;
                    o[2][d] = nA ? MA[d] / nA : 0; o[3][d] = nB ? MB[d] / nB : 0;
                }
                for (int q = 0; q < 4; q++) { double m2 = o[q][0]*o[q][0] + o[q][1]*o[q][1] + o[q][2]*o[q][2]; acc[2*q] += m2; acc[2*q+1] += m2*m2; acc[8+q] += sqrt(m2); }
                if (haspos) {
                    double sk = 0;
                    for (int d = 0; d < 3; d++) {
                        double re[3] = {0}, im[3] = {0};
                        for (int i = 0; i < N; i++) for (int c = 0; c < 3; c++) {
                            re[c] += e[3*i+c] * cs[6*i+2*d]; im[c] += e[3*i+c] * cs[6*i+2*d+1]; }
                        for (int c = 0; c < 3; c++) sk += re[c]*re[c] + im[c]*im[c];
                    }
                    skacc += sk / (3.0 * N * (double)N);
                }
                eacc += en / N; nm++;
            }
        }
        printf("%.6f %.8f", T, eacc / nm);
        for (int q = 0; q < 12; q++) printf(" %.10e", acc[q] / nm);
        if (haspos) printf(" %.10e", skacc / nm);
        printf("\n"); fflush(stdout);
    }
    return 0;
}
