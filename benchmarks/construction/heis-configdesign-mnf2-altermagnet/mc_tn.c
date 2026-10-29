/* Classical Heisenberg Monte Carlo for the two-sublattice MnF2 exchange model (reference T_N).
 *
 *   gcc -O3 -o mc_tn mc_tn.c -lm
 *   ./mc_tn bonds.txt L T_meV n_equil n_meas seed > out.txt
 *
 * H = - sum_{i != j} J_ij e_i.e_j over ordered pairs (unit vectors), so the local field is
 * h_i = 2 sum_j J_ij e_j. Lattice: L x L x L conventional rutile cells, 2 Mn per cell (s = 0 at (1/2,1/2,1/2),
 * s = 1 at the origin), periodic. bonds.txt: lines "s ds dx dy dz J" (neighbour of sublattice s in sublattice
 * ds, in the cell displaced by (dx,dy,dz)). One sweep = one heat-bath sweep + 2 over-relaxation sweeps
 * (random site order). Output per measurement block: <m_st^2>, <m_st^4>, <E>, <E^2> with m_st the staggered
 * magnetisation per site and E the energy per Mn (meV).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned long long rng_s[2];
static inline unsigned long long rotl(const unsigned long long x, int k) { return (x << k) | (x >> (64 - k)); }
static inline unsigned long long next64(void) {
    unsigned long long s0 = rng_s[0], s1 = rng_s[1], r = s0 + s1;
    s1 ^= s0; rng_s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); rng_s[1] = rotl(s1, 36);
    return r;
}
static inline double urand(void) { return (next64() >> 11) * (1.0 / 9007199254740992.0); }

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: mc_tn bonds L T n_equil n_meas seed\n"); return 1; }
    int L = atoi(argv[2]); double T = atof(argv[3]); long neq = atol(argv[4]), nmeas = atol(argv[5]);
    rng_s[0] = 0x9E3779B97F4A7C15ULL ^ (unsigned long long)atoll(argv[6]); rng_s[1] = 0xD1B54A32D192ED03ULL;
    for (int i = 0; i < 20; i++) next64();
    FILE *f = fopen(argv[1], "r");
    int bs[512], bds[512], bx[512], by[512], bz[512]; double bJ[512]; int nb = 0;
    while (nb < 512 && fscanf(f, "%d %d %d %d %d %lf", &bs[nb], &bds[nb], &bx[nb], &by[nb], &bz[nb], &bJ[nb]) == 6) nb++;
    fclose(f);
    int N = 2 * L * L * L;
    int maxnb = 0, cnt[2] = {0, 0};
    for (int b = 0; b < nb; b++) cnt[bs[b]]++;
    maxnb = cnt[0] > cnt[1] ? cnt[0] : cnt[1];
    int *nbr = malloc(sizeof(int) * N * maxnb); double *nJ = malloc(sizeof(double) * N * maxnb); int *nn = calloc(N, sizeof(int));
    for (int x = 0; x < L; x++) for (int y = 0; y < L; y++) for (int z = 0; z < L; z++) for (int s = 0; s < 2; s++) {
        int i = ((x * L + y) * L + z) * 2 + s;
        for (int b = 0; b < nb; b++) if (bs[b] == s) {
            int X = ((x + bx[b]) % L + L) % L, Y = ((y + by[b]) % L + L) % L, Z = ((z + bz[b]) % L + L) % L;
            int j = ((X * L + Y) * L + Z) * 2 + bds[b];
            nbr[i * maxnb + nn[i]] = j; nJ[i * maxnb + nn[i]] = bJ[b]; nn[i]++;
        }
    }
    double *e = malloc(sizeof(double) * 3 * N);
    for (int i = 0; i < N; i++) { double sg = (i % 2) ? -1.0 : 1.0; e[3*i] = 0; e[3*i+1] = 0; e[3*i+2] = sg; }
    double beta = 1.0 / T;
    long total = neq + nmeas; int block = 100;
    double am2 = 0, am4 = 0, aE = 0, aE2 = 0; long nacc = 0;
    for (long sweep = 0; sweep < total; sweep++) {
        for (int pass = 0; pass < 3; pass++) {
            for (int k = 0; k < N; k++) {
                int i = (int)(urand() * N);
                double hx = 0, hy = 0, hz = 0;
                for (int b = 0; b < nn[i]; b++) { int j = nbr[i * maxnb + b]; double J = nJ[i * maxnb + b];
                    hx += J * e[3*j]; hy += J * e[3*j+1]; hz += J * e[3*j+2]; }
                hx *= 2; hy *= 2; hz *= 2;
                double h = sqrt(hx*hx + hy*hy + hz*hz);
                if (h < 1e-12) continue;
                double ux = hx / h, uy = hy / h, uz = hz / h;
                if (pass == 0) {   /* heat bath: p(c) ~ exp(beta h c) */
                    double a = beta * h, c;
                    double r = urand();
                    if (a > 1e-8) c = 1.0 + log(r + (1.0 - r) * exp(-2.0 * a)) / a; else c = 2.0 * r - 1.0;
                    if (c > 1) c = 1; if (c < -1) c = -1;
                    double sn = sqrt(1.0 - c * c), ph = 2 * M_PI * urand();
                    /* orthonormal frame around u */
                    double vx, vy, vz;
                    if (fabs(uz) < 0.9) { vx = -uy; vy = ux; vz = 0; } else { vx = 0; vy = -uz; vz = uy; }
                    double vn = sqrt(vx*vx + vy*vy + vz*vz); vx /= vn; vy /= vn; vz /= vn;
                    double wx = uy*vz - uz*vy, wy = uz*vx - ux*vz, wz = ux*vy - uy*vx;
                    e[3*i] = c*ux + sn*(cos(ph)*vx + sin(ph)*wx);
                    e[3*i+1] = c*uy + sn*(cos(ph)*vy + sin(ph)*wy);
                    e[3*i+2] = c*uz + sn*(cos(ph)*vz + sin(ph)*wz);
                } else {           /* over-relaxation */
                    double d = e[3*i]*ux + e[3*i+1]*uy + e[3*i+2]*uz;
                    e[3*i] = 2*d*ux - e[3*i]; e[3*i+1] = 2*d*uy - e[3*i+1]; e[3*i+2] = 2*d*uz - e[3*i+2];
                }
            }
        }
        if (sweep >= neq) {
            double mx = 0, my = 0, mz = 0, E = 0;
            for (int i = 0; i < N; i++) { double sg = (i % 2) ? -1.0 : 1.0; mx += sg*e[3*i]; my += sg*e[3*i+1]; mz += sg*e[3*i+2];
                for (int b = 0; b < nn[i]; b++) { int j = nbr[i * maxnb + b];
                    E -= nJ[i * maxnb + b] * (e[3*i]*e[3*j] + e[3*i+1]*e[3*j+1] + e[3*i+2]*e[3*j+2]); } }
            double m2 = (mx*mx + my*my + mz*mz) / ((double)N * N); E /= N;
            am2 += m2; am4 += m2 * m2; aE += E; aE2 += E * E; nacc++;
            if (nacc == block) { printf("%.10g %.10g %.10g %.10g\n", am2 / nacc, am4 / nacc, aE / nacc, aE2 / nacc); fflush(stdout);
                am2 = am4 = aE = aE2 = 0; nacc = 0; }
        }
    }
    return 0;
}
