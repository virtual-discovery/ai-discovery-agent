/* Rejection-free vacancy kinetic Monte Carlo on a periodic bcc lattice (single vacancy).
 *
 * Barrier of a hop (vacancy at v exchanges with atom X at n):
 *     Eb = KRA(X, environment of the pair) + (E_after - E_before)/2
 * KRA: K0[X] + sum_orbits sum_sites W[X][o][species(site)] / n_sites(o)   (orbit = class of sites around the pair)
 *      + sum_ij Q[X][i][j] c_i c_j, c = composition of the first n_qorb orbits (the nearest-neighbour shell of the pair)
 * End states: pair Hamiltonian E = sum_shells sum_pairs V[s][a][b] (species 0-5 atoms, 6 vacancy).
 * Optional canonical (Kawasaki, non-local swap) Monte Carlo equilibration with the same pair Hamiltonian.
 *
 * Reads a whitespace-separated input from stdin (see kmc_io.py) and writes one line per trajectory.
 * Positions are integers in units of a/2; times in units of 1/nu.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#define NSP 6
#define VAC 6
#define MAXORB 16
#define MAXSH 6
#define MAXCK 64

static uint64_t rs[4];
static inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t next64(void) {
    const uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45);
    return r;
}
static double urand(void) { return ((next64() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
static void seed_rng(uint64_t s) {
    for (int i = 0; i < 4; i++) { s += 0x9E3779B97F4A7C15ULL; uint64_t z = s;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; rs[i] = z ^ (z >> 31); }
}

static int L, N, M2;
static int site_index(int x, int y, int z) {
    x = ((x % M2) + M2) % M2; y = ((y % M2) + M2) % M2; z = ((z % M2) + M2) % M2;
    int b = x & 1; int i = (x - b) / 2, j = (y - b) / 2, k = (z - b) / 2;
    return ((i * L + j) * L + k) * 2 + b;
}
static int pos[1 << 20][3];

static int nsteps, ntraj, nmc, nck, mode_thermo, n_orb, n_sh;
static long seed;
static double T, nu, kT;
static int counts[NSP];
static double K0[NSP], W[NSP][MAXORB][NSP + 1], V[MAXSH][NSP + 1][NSP + 1], Q[NSP][NSP][NSP];
static int n_qorb;
static int orb_n[MAXORB], *orb_off[MAXORB];
static int sh_n[MAXSH], *sh_off[MAXSH];
static int *orbtab;   /* [N][8][tot_orb_sites] */
static int tot_orb;
static int *shtab;    /* [N][tot_sh] */
static int tot_sh;
static const int dirs[8][3] = {{-1,-1,-1},{-1,-1,1},{-1,1,-1},{-1,1,1},{1,-1,-1},{1,-1,1},{1,1,-1},{1,1,1}};
static int nntab[1 << 20][8];

static unsigned char *occ;
static int *disp;   /* displacement (a/2 units) of the occupant of each site, 3 per site */

static double site_e(int i) {
    double e = 0; int a = occ[i]; int *nb = shtab + (size_t)i * tot_sh; int c = 0;
    for (int s = 0; s < n_sh; s++) for (int q = 0; q < sh_n[s]; q++, c++) e += V[s][a][occ[nb[c]]];
    return e;
}
static double bond_e(int i, int j) {
    double e = 0; int *nb = shtab + (size_t)i * tot_sh; int c = 0;
    for (int s = 0; s < n_sh; s++) for (int q = 0; q < sh_n[s]; q++, c++) if (nb[c] == j) e += V[s][occ[i]][occ[j]];
    return e;
}
/* energy change for swapping the occupants of sites i and j (pair Hamiltonian) */
static double swap_de(int i, int j) {
    if (mode_thermo == 0) return 0.0;
    double before = site_e(i) + site_e(j) - bond_e(i, j);
    unsigned char t = occ[i]; occ[i] = occ[j]; occ[j] = t;
    double after = site_e(i) + site_e(j) - bond_e(i, j);
    t = occ[i]; occ[i] = occ[j]; occ[j] = t;
    return after - before;
}
static double kra(int v, int d, int X) {
    double e = K0[X]; int *st = orbtab + ((size_t)v * 8 + d) * tot_orb; int c = 0;
    double nn[NSP + 1] = {0}; int nq = 0;
    for (int o = 0; o < n_orb; o++) { double acc = 0;
        for (int q = 0; q < orb_n[o]; q++, c++) { int sg = occ[st[c]]; acc += W[X][o][sg]; if (o < n_qorb) { nn[sg] += 1; nq++; } }
        e += acc / orb_n[o]; }
    if (n_qorb > 0) for (int i = 0; i < NSP; i++) for (int j = 0; j < NSP; j++) e += Q[X][i][j] * nn[i] * nn[j] / ((double)nq * nq);
    return e;
}

int main(void) {
#define RD(fmt, p) if (scanf(fmt, p) != 1) { fprintf(stderr, "bad input\n"); return 2; }
    RD("%d", &L); RD("%d", &nsteps); RD("%d", &ntraj); RD("%ld", &seed); RD("%lf", &T); RD("%lf", &nu);
    RD("%d", &nmc); RD("%d", &nck); RD("%d", &mode_thermo);
    for (int s = 0; s < NSP; s++) RD("%d", &counts[s]);
    for (int s = 0; s < NSP; s++) RD("%lf", &K0[s]);
    RD("%d", &n_orb);
    if (n_orb > MAXORB || nck > MAXCK) { fprintf(stderr, "too many orbits/checkpoints\n"); return 2; }
    tot_orb = 0;
    for (int o = 0; o < n_orb; o++) { RD("%d", &orb_n[o]); orb_off[o] = malloc(sizeof(int) * 3 * orb_n[o]);
        for (int q = 0; q < 3 * orb_n[o]; q++) RD("%d", &orb_off[o][q]); tot_orb += orb_n[o]; }
    for (int s = 0; s < NSP; s++) for (int o = 0; o < n_orb; o++) for (int e = 0; e < NSP; e++) RD("%lf", &W[s][o][e]);
    RD("%d", &n_qorb);
    if (n_qorb > n_orb) { fprintf(stderr, "quadratic orbits exceed orbits\n"); return 2; }
    for (int s = 0; s < NSP; s++) for (int i = 0; i < NSP; i++) for (int j = 0; j < NSP; j++) RD("%lf", &Q[s][i][j]);
    RD("%d", &n_sh);
    if (n_sh > MAXSH) { fprintf(stderr, "too many shells\n"); return 2; }
    tot_sh = 0;
    for (int s = 0; s < n_sh; s++) { RD("%d", &sh_n[s]); sh_off[s] = malloc(sizeof(int) * 3 * sh_n[s]);
        for (int q = 0; q < 3 * sh_n[s]; q++) RD("%d", &sh_off[s][q]); tot_sh += sh_n[s]; }
    for (int s = 0; s < n_sh; s++) for (int a = 0; a <= NSP; a++) for (int b = 0; b <= NSP; b++) RD("%lf", &V[s][a][b]);
    M2 = 2 * L; N = 2 * L * L * L; kT = 8.617333262e-5 * T;
    int nat = 0; for (int s = 0; s < NSP; s++) nat += counts[s];
    if (nat != N - 1 || N > (1 << 20)) { fprintf(stderr, "composition must fill N-1 sites\n"); return 2; }
    for (int i = 0; i < L; i++) for (int j = 0; j < L; j++) for (int k = 0; k < L; k++) for (int b = 0; b < 2; b++) {
        int n = ((i * L + j) * L + k) * 2 + b; pos[n][0] = 2 * i + b; pos[n][1] = 2 * j + b; pos[n][2] = 2 * k + b; }
    orbtab = malloc(sizeof(int) * (size_t)N * 8 * (tot_orb ? tot_orb : 1));
    shtab = malloc(sizeof(int) * (size_t)N * (tot_sh ? tot_sh : 1));
    for (int n = 0; n < N; n++) {
        for (int d = 0; d < 8; d++) {
            nntab[n][d] = site_index(pos[n][0] + dirs[d][0], pos[n][1] + dirs[d][1], pos[n][2] + dirs[d][2]);
            int c = 0;
            for (int o = 0; o < n_orb; o++) for (int q = 0; q < orb_n[o]; q++, c++) {
                int *f = orb_off[o] + 3 * q;   /* offsets for direction (+,+,+); flip signs per axis */
                orbtab[((size_t)n * 8 + d) * tot_orb + c] = site_index(pos[n][0] + dirs[d][0] * f[0],
                    pos[n][1] + dirs[d][1] * f[1], pos[n][2] + dirs[d][2] * f[2]);
            }
        }
        int c = 0;
        for (int s = 0; s < n_sh; s++) for (int q = 0; q < sh_n[s]; q++, c++) {
            int *f = sh_off[s] + 3 * q;
            shtab[(size_t)n * tot_sh + c] = site_index(pos[n][0] + f[0], pos[n][1] + f[1], pos[n][2] + f[2]);
        }
    }
    occ = malloc(N); disp = malloc(sizeof(int) * 3 * (size_t)N);
    seed_rng((uint64_t)seed);
    for (int tr = 0; tr < ntraj; tr++) {
        /* random alloy with exact counts, vacancy at a random site */
        int c = 0; for (int s = 0; s < NSP; s++) for (int q = 0; q < counts[s]; q++) occ[c++] = s;
        occ[c] = VAC;
        for (int i = N - 1; i > 0; i--) { int j = (int)(urand() * (i + 1)); unsigned char t = occ[i]; occ[i] = occ[j]; occ[j] = t; }
        /* canonical Monte Carlo: non-local swaps of any two different occupants */
        long acc_mc = 0;
        for (long it = 0; it < (long)nmc * N; it++) {
            int i = (int)(urand() * N), j = (int)(urand() * N);
            if (occ[i] == occ[j]) continue;
            double de = swap_de(i, j);
            if (de <= 0 || urand() < exp(-de / kT)) { unsigned char t = occ[i]; occ[i] = occ[j]; occ[j] = t; acc_mc++; }
        }
        int v = -1; for (int i = 0; i < N; i++) if (occ[i] == VAC) v = i;
        { long pc[NSP][NSP] = {{0}};   /* nearest-neighbour pair counts after equilibration (ordered pairs) */
          for (int i = 0; i < N; i++) if (occ[i] != VAC) for (int d = 0; d < 8; d++) { int j = nntab[i][d]; if (occ[j] != VAC) pc[occ[i]][occ[j]]++; }
          printf("S %d", tr); for (int a = 0; a < NSP; a++) for (int b = 0; b < NSP; b++) printf(" %ld", pc[a][b]); printf("\n"); }
        memset(disp, 0, sizeof(int) * 3 * (size_t)N);
        double t = 0, sum_enc = 0, sum_acc = 0; long hops[NSP] = {0};
        int vd[3] = {0, 0, 0};
        int ck = 0;
        for (int step = 1; step <= nsteps; step++) {
            double rate[8], eb[8], tot = 0;
            for (int d = 0; d < 8; d++) {
                int n = nntab[v][d]; int X = occ[n];
                double b = kra(v, d, X) + 0.5 * swap_de(v, n);
                eb[d] = b; rate[d] = exp(-b / kT); tot += rate[d]; sum_enc += b;
            }
            double r = urand() * tot, cum = 0; int d = 0;
            for (d = 0; d < 7; d++) { cum += rate[d]; if (r < cum) break; }
            int n = nntab[v][d]; int X = occ[n];
            sum_acc += eb[d]; hops[X]++;
            t += -log(urand()) / tot;
            /* atom moves n -> v (displacement -dirs[d]), vacancy moves v -> n */
            disp[3 * v + 0] = disp[3 * n + 0] - dirs[d][0];
            disp[3 * v + 1] = disp[3 * n + 1] - dirs[d][1];
            disp[3 * v + 2] = disp[3 * n + 2] - dirs[d][2];
            disp[3 * n + 0] = disp[3 * n + 1] = disp[3 * n + 2] = 0;
            occ[v] = X; occ[n] = VAC;
            vd[0] += dirs[d][0]; vd[1] += dirs[d][1]; vd[2] += dirs[d][2];
            v = n;
            if (nck > 0 && step == (long)nsteps * (ck + 1) / nck) {
                double r2[NSP] = {0};
                for (int i = 0; i < N; i++) if (occ[i] != VAC) {
                    int *q = disp + 3 * i; r2[occ[i]] += (double)q[0] * q[0] + (double)q[1] * q[1] + (double)q[2] * q[2]; }
                printf("C %d %d %d %.10e %.10e", tr, ck, step, t / nu, (double)vd[0] * vd[0] + (double)vd[1] * vd[1] + (double)vd[2] * vd[2]);
                for (int s = 0; s < NSP; s++) printf(" %.10e", r2[s]);
                printf("\n");
                ck++;
            }
        }
        double r2[NSP] = {0};
        for (int i = 0; i < N; i++) if (occ[i] != VAC) {
            int *q = disp + 3 * i; r2[occ[i]] += (double)q[0] * q[0] + (double)q[1] * q[1] + (double)q[2] * q[2]; }
        printf("T %d %.10e %.10e", tr, t / nu, (double)vd[0] * vd[0] + (double)vd[1] * vd[1] + (double)vd[2] * vd[2]);
        for (int s = 0; s < NSP; s++) printf(" %.10e", r2[s]);
        for (int s = 0; s < NSP; s++) printf(" %ld", hops[s]);
        printf(" %.8f %.8f %.6f\n", sum_enc / (8.0 * nsteps), sum_acc / nsteps, nmc ? (double)acc_mc / ((double)nmc * N) : 0.0);
        fflush(stdout);
    }
    return 0;
}
